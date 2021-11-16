#ifndef TEMPERATURE_H
#define TEMPERATURE_H

// 每個生態域均有一個溫度數值來決定該位置是否能夠降水。該數值低於0.15即下雪，0.15-0.95即下雨，高於1.0將會使區域保持乾旱。這些數值也用於決定不同生態域中積雪的高度。每從預設海平面（Y=64）向上升高1米，溫度將會降低0.0016（1⁄625）。海平面以下不會有溫度變化。

class Temperature
{
public:
    Temperature();
    ~Temperature();

    Temperature *getInstance() {
        if (instance == nullptr) {
            instance = new Temperature();
        }
        return instance;
    }
    static float getTemperature(int level) {
        if (level <= m_nSeaLevel) return 0;
        return 0 - (level - m_nSeaLevel) * m_fTempStep;
    }

private:
    static Temperature *instance;
    static int m_nSeaLevel;
    static const float m_fTempStep;
};

#endif // TEMPERATURE_H
