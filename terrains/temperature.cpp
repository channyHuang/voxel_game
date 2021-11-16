#include "temperature.h"

Temperature* Temperature::instance = nullptr;
int Temperature::m_nSeaLevel = 0;
const float Temperature::m_fTempStep = 1.f / 625.f;

Temperature::Temperature()
{

}
