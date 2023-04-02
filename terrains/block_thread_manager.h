#ifndef VOXEL_BLOCK_THREAD_MANAGER_H
#define VOXEL_BLOCK_THREAD_MANAGER_H

#include <QDebug>

#include <algorithm>
#include <functional>
#include <vector>
#include <mutex>
#include <thread>
//#include <semaphore>
#include <unordered_map>

#include "boxi.h"
#include "vector3i.h"
#include "utility.h"
#include "array_slice.h"
#include "fixed_array.h"
#include "counting_time.h"

typedef Boxi Rect3i;

template <typename InputBlockData_T, typename OutputBlockData_T>
class VoxelBlockThreadManager {
public:
    static const int MAX_JOBS = 8;

    struct InputBlock {
        InputBlockData_T data;
        Vector3i position;
        uint8_t lod = 0;
        bool can_be_discarded = true;
        float sort_heuristic = 0;
    };

    struct OutputBlock {
        OutputBlockData_T data;
        Vector3i position;
        uint8_t lod = 0;
        bool drop_hint = false;
    };

    struct Input {
        std::vector<InputBlock> blocks;
        Vector3i priority_position;
        Vector3 priority_direction;
        int exclusive_region_extent = 0;
        int exclusive_region_max_lod = Math::MAX_LOD;
        bool use_exclusive_region = false;
        int max_lod_index = 0;

        bool is_empty() const {
            return blocks.empty();
        }
    };

    struct ProcessorStats {
        int file_openings = 0;
        int time_spent_opening_files = 0;
    };

    struct Stats {
        bool first = true;
        uint64_t min_time = 0;
        uint64_t max_time = 0;
        uint64_t sorting_time = 0;
        uint32_t remaining_blocks[MAX_JOBS];
        uint32_t thread_count = 0;
        uint32_t dropped_count = 0;
        ProcessorStats processor;

        Stats() {
            for (int i = 0; i < MAX_JOBS; ++i) {
                remaining_blocks[i] = 0;
            }
        }
    };

    struct Output {
        std::vector<OutputBlock> blocks;
        Stats stats;
    };

    typedef std::function<void(ArraySlice<InputBlock>, ArraySlice<OutputBlock>, ProcessorStats &)> BlockProcessingFunc;

    VoxelBlockThreadManager(
            unsigned int job_count,
            unsigned int sync_interval_ms,
            ArraySlice<BlockProcessingFunc> processors,
            bool duplicate_rejection = true,
            unsigned int batch_count = 1) {

        _job_count = job_count;

        for (unsigned int i = 0; i < MAX_JOBS; ++i) {
            JobData &job = _jobs[i];
            job.job_index = i;
            job.duplicate_rejection = duplicate_rejection;
            job.sync_interval_ms = sync_interval_ms;
            job.batch_count = batch_count;
        }

        for (unsigned int i = 0; i < _job_count; ++i) {
            JobData &job = _jobs[i];

            //job.input_mutex.unlock();
            //job.output_mutex.unlock();
            //job.sema.release();
            job.thread = std::thread(_thread_func, &job);
            job.needs_sort = true;
            job.processor = processors[i];
        }
    }

    ~VoxelBlockThreadManager() {
        for (unsigned int i = 0; i < _job_count; ++i) {
            JobData &job = _jobs[i];
            job.thread_exit = true;
            //job.sema.release();
        }

        for (unsigned int i = 0; i < _job_count; ++i) {
            JobData &job = _jobs[i];

            job.thread.join();

            //job.sema.release();
            //job.input_mutex.unlock();
            //job.output_mutex.unlock();
        }
    }

    void push(const Input &input) {
        unsigned int replaced_blocks = 0;
        unsigned int highest_pending_count = 0;
        unsigned int lowest_pending_count = 0;

        for (unsigned int job_index = 0; job_index < _job_count; ++job_index) {
            JobData &job = _jobs[job_index];

            job.input_mutex.lock();

            highest_pending_count = std::max(highest_pending_count, (unsigned int)job.shared_input.blocks.size());
            lowest_pending_count = std::min(lowest_pending_count, (unsigned int)job.shared_input.blocks.size());
        }

        unsigned int i = 0;

        unsigned int median_pending_count = lowest_pending_count + (highest_pending_count - lowest_pending_count) / 2;

        for (unsigned int job_index = 0; job_index < _job_count && i < input.blocks.size(); ++job_index) {
            JobData &job = _jobs[job_index];
            unsigned int pending_count = job.shared_input.blocks.size();

            unsigned int count = std::min(median_pending_count - pending_count, (unsigned int)input.blocks.size());

            if (count > 0) {
                if (i + count > input.blocks.size()) {
                    count = input.blocks.size() - i;
                }
                replaced_blocks += push_block_requests(job, input.blocks, i, count);
                i += count;
            }
        }

        unsigned int base_count = (input.blocks.size() - i) / _job_count;
        unsigned int remainder = (input.blocks.size() - i) % _job_count;
        for (unsigned int job_index = 0; job_index < _job_count && i < input.blocks.size(); ++job_index) {
            JobData &job = _jobs[job_index];

            unsigned int count = base_count;
            if (remainder > 0) {
                ++count;
                --remainder;
            }

            if (i + count > input.blocks.size()) {
                replaced_blocks += push_block_requests(job, input.blocks, i, input.blocks.size() - i);
            } else {
                replaced_blocks += push_block_requests(job, input.blocks, i, count);
                i += count;
            }
        }

        for (unsigned int job_index = 0; job_index < _job_count; ++job_index) {
            JobData &job = _jobs[job_index];

            if (job.shared_input.priority_position != input.priority_position || input.blocks.size() > 0) {
                job.needs_sort = true;
            }

            job.shared_input.priority_position = input.priority_position;

            if (input.use_exclusive_region) {
                job.shared_input.use_exclusive_region = true;
                job.shared_input.exclusive_region_extent = input.exclusive_region_extent;
                job.shared_input.exclusive_region_max_lod = input.exclusive_region_max_lod;
            }

            bool should_run = !job.shared_input.is_empty();

            job.input_mutex.unlock();

            if (should_run) {
                 //job.sema.release();
            }
        }

        if (replaced_blocks > 0) {
            char str[1024];
            sprintf(str, "VoxelBlockProcessor: {1} blocks already in queue were replaced");
        }

    }

    void pop(Output &output) {
        output.stats = Stats();
        output.stats.thread_count = _job_count;

        for (unsigned int i = 0; i < _job_count; ++i) {

            JobData &job = _jobs[i];
            {
                std::lock_guard<std::mutex>(job.output_mutex);

                output.blocks.insert(output.blocks.end(), job.shared_output.blocks.begin(), job.shared_output.blocks.end());
                merge_stats(output.stats, job.shared_output.stats, i);
                job.shared_output.blocks.clear();
            }
        }
    }

private:
    struct JobData {
        Input shared_input;
        Output shared_output;

        std::mutex input_mutex;
        std::mutex output_mutex;

        FixedArray<std::unordered_map<Vector3i, int, Vector3iHasher>, Math::MAX_LOD> shared_input_block_indexes;
        bool needs_sort = false;

        bool thread_exit = false;
        //------------------------

        Input input;
        Output output;

        //std::binary_semaphore sema = std::binary_semaphore(1);
        std::thread thread;
        uint32_t sync_interval_ms = 100;
        uint32_t job_index = -1;
        bool duplicate_rejection = false;
        int batch_count = 0;

        BlockProcessingFunc processor;
    };

    static void merge_stats(Stats &a, const Stats &b, int job_index) {
        a.max_time = std::max(a.max_time, b.max_time);
        a.min_time = std::min(a.min_time, b.min_time);
        a.remaining_blocks[job_index] = b.remaining_blocks[job_index];
        a.sorting_time += b.sorting_time;
        a.dropped_count += b.dropped_count;

        a.processor.file_openings += b.processor.file_openings;
        a.processor.time_spent_opening_files += b.processor.time_spent_opening_files;
    }

    unsigned int push_block_requests(JobData &job, const std::vector<InputBlock> &input_blocks, int begin, int count) {
        unsigned int replaced_blocks = 0;
        unsigned int end = begin + count;

        for (unsigned int i = begin; i < end; ++i) {
            const InputBlock &block = input_blocks[i];

            if (job.duplicate_rejection) {
                int index = job.shared_input_block_indexes[block.lod][block.position];

                if (index) {
                    ++replaced_blocks;
                    job.shared_input.blocks[index] = block;

                } else {
                    unsigned int j = job.shared_input.blocks.size();
                    job.shared_input.blocks.push_back(block);
                    job.shared_input_block_indexes[block.lod][block.position] = j;
                }

            } else {
                job.shared_input.blocks.push_back(block);
            }
        }

        return replaced_blocks;
    }

    static void _thread_func(void *p_data) {
        JobData *data = reinterpret_cast<JobData *>(p_data);
        thread_func(*data);
    }

    static void thread_func(JobData &data) {
        while (true) {
            uint32_t sync_time = CountingTime::getInstance()->get_ticks_msec() + data.sync_interval_ms;

            unsigned int queue_index = 0;
            Stats stats;

            bool should_exit = false;
            if (data.thread_exit) {
                should_exit = true;
            }

            thread_sync(data, queue_index, stats, stats.sorting_time, stats.dropped_count);

            while (!data.input.blocks.empty()) {
                if (!data.input.blocks.empty()) {
                    if (should_exit) {
                        shift_up(data.input.blocks, queue_index);
                        queue_index = 0;

                        unordered_remove_if(data.input.blocks,
                                [](const InputBlock &b) {
                                    return b.can_be_discarded;
                                });
                    }

                    unsigned int input_begin = queue_index;
                    unsigned int batch_count = data.batch_count;

                    if (input_begin + batch_count > data.input.blocks.size()) {
                        batch_count = data.input.blocks.size() - input_begin;
                    }

                    if (batch_count > 0) {
                        uint64_t time_before = CountingTime::getInstance()->get_ticks_usec();

                        unsigned int output_begin = data.output.blocks.size();
                        data.output.blocks.resize(data.output.blocks.size() + batch_count);

                        for (unsigned int i = 0; i < batch_count; ++i) {
                            InputBlock &ib = data.input.blocks[input_begin + i];
                            OutputBlock &ob = data.output.blocks[output_begin + i];
                            ob.position = ib.position;
                            ob.lod = ib.lod;
                        }

                        data.processor(
                                ArraySlice<InputBlock>(data.input.blocks, input_begin, input_begin + batch_count),
                                ArraySlice<OutputBlock>(&data.output.blocks[0], output_begin, output_begin + batch_count),
                                stats.processor);

                        uint64_t time_taken = (CountingTime::getInstance()->get_ticks_usec() - time_before) / batch_count;

                        // Do some stats
                        if (stats.first) {
                            stats.first = false;
                            stats.min_time = time_taken;
                            stats.max_time = time_taken;
                        } else {
                            if (time_taken < stats.min_time) {
                                stats.min_time = time_taken;
                            }
                            if (time_taken > stats.max_time) {
                                stats.max_time = time_taken;
                            }
                        }
                    }

                    queue_index += batch_count;
                    if (queue_index >= data.input.blocks.size()) {
                        data.input.blocks.clear();
                    }
                }

                uint32_t time = CountingTime::getInstance()->get_ticks_msec();
                if (time >= sync_time || data.input.blocks.empty()) {
                    uint64_t sort_time;
                    unsigned int dropped_count;
                    thread_sync(data, queue_index, stats, sort_time, dropped_count);

                    sync_time = CountingTime::getInstance()->get_ticks_msec() + data.sync_interval_ms;
                    queue_index = 0;
                    stats = Stats();
                    stats.sorting_time = sort_time;
                    stats.dropped_count = dropped_count;
                }
            }

            if (should_exit) {
                break;
            }

            // Wait for future wake-up
            //data.sema.acquire();
        }

        qDebug() << "exit";
    }

    static inline float get_priority_heuristic(const InputBlock &a, const Vector3i &viewer_block_pos, const Vector3 &viewer_direction, int max_lod) {
        int f = 1 << a.lod;
        Vector3i p = a.position * f;
        float d = std::sqrt(p.distance(viewer_block_pos) + 0.1f);
        float dp = viewer_direction.dot(viewer_block_pos.to_vec3() / d);
        return (max_lod - a.lod) * 10000.f + d + (1.f - dp) * 4.f * f;
    }

    struct BlockUpdateComparator {
        inline bool operator()(const InputBlock &a, const InputBlock &b) const {
            return a.sort_heuristic < b.sort_heuristic;
        }
    };

    static void thread_sync(JobData &data, unsigned int queue_index, Stats stats, uint64_t &out_sort_time, unsigned int &out_dropped_count) {
        if (!data.input.blocks.empty()) {

            if (queue_index >= data.input.blocks.size()) {
                data.input.blocks.clear();

            } else if (queue_index > 0) {

                shift_up(data.input.blocks, queue_index);
            }
        }

        stats.remaining_blocks[data.job_index] = data.input.blocks.size();
        bool needs_sort;

        {
            std::lock_guard<std::mutex>(data.input_mutex);

            append_array(data.input.blocks, data.shared_input.blocks);

            data.input.priority_position = data.shared_input.priority_position;

            if (data.shared_input.use_exclusive_region) {
                data.input.use_exclusive_region = true;
                data.input.exclusive_region_extent = data.shared_input.exclusive_region_extent;
                data.input.exclusive_region_max_lod = data.shared_input.exclusive_region_max_lod;
            }

            data.shared_input.blocks.clear();

            if (data.duplicate_rejection) {
                for (unsigned int lod_index = 0; lod_index < data.shared_input_block_indexes.size(); ++lod_index) {
                    data.shared_input_block_indexes[lod_index].clear();
                }
            }

            needs_sort = data.needs_sort;
            data.needs_sort = false;
        }

        if (!data.output.blocks.empty()) {
            std::lock_guard<std::mutex>(data.output_mutex);

            data.shared_output.blocks.insert(data.shared_output.blocks.end(), data.output.blocks.begin(), data.output.blocks.end());
            data.shared_output.stats = stats;
            data.output.blocks.clear();
        }

        int dropped_count = 0;
        if (data.input.use_exclusive_region) {
            for (unsigned int i = 0; i < data.input.blocks.size(); ++i) {
                const InputBlock &ib = data.input.blocks[i];

                if (!ib.can_be_discarded || ib.lod >= data.input.exclusive_region_max_lod) {
                    continue;
                }

                Rect3i box = Rect3i::from_center_extents(data.input.priority_position >> ib.lod, Vector3i(data.input.exclusive_region_extent));

                if (!box.contains(ib.position)) {
                    OutputBlock ob;
                    ob.position = ib.position;
                    ob.lod = ib.lod;
                    ob.drop_hint = true;
                    data.output.blocks.push_back(ob);

                    const InputBlock &shifted_block = data.input.blocks.back();

                    data.input.blocks[i] = shifted_block;
                    data.input.blocks.pop_back();

                    --i;

                    ++dropped_count;
                }
            }
        }

        if (dropped_count > 0) {
            char str[1024];
            sprintf(str, "Dropped {0} blocks from thread");
            out_dropped_count = dropped_count;
        }

        uint64_t time_before = CountingTime::getInstance()->get_ticks_usec();

        if (!data.input.blocks.empty() && needs_sort) {
            for (auto it = data.input.blocks.begin(); it != data.input.blocks.end(); ++it) {
                InputBlock &ib = *it;
                ib.sort_heuristic = get_priority_heuristic(ib,
                        data.input.priority_position,
                        data.input.priority_direction,
                        data.input.max_lod_index);
            }

            auto BlockUpdateComparator = [](const InputBlock &a, const InputBlock &b) {
                return a.sort_heuristic < b.sort_heuristic;
            };

            std::sort(data.input.blocks.begin(), data.input.blocks.end(), BlockUpdateComparator);
        }

        out_sort_time = CountingTime::getInstance()->get_ticks_usec() - time_before;
    }

    JobData _jobs[MAX_JOBS];
    unsigned int _job_count = 0;
};

#endif // VOXEL_BLOCK_THREAD_MANAGER_H
