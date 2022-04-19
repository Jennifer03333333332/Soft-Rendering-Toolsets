
#include "../geometry/spline.h"
#include "debug.h"

template<typename T>
T Spline<T>::cubic_unit_spline(float time, const T& position0, const T& position1,
                               const T& tangent0, const T& tangent1) {

    // TODO (Animation): Task 1a
    // Given time in [0,1] compute the cubic spline coefficients and use them to compute
    // the interpolated value at time 'time' based on the positions & tangents
    float t_squared = time*time;
    float t_cubic = t_squared*time;

    float h00 = 2.0f*t_cubic - 3.0f*t_squared + 1.0f;
    float h10 = t_cubic - 2.0f*t_squared + time;
    float h01 = -2.0f*t_cubic + 3.0f*t_squared;
    float h11 = t_cubic - t_squared;
    return h00*position0 + h10*tangent0 + h01*position1 + h11*tangent1;

    //return T();
}

template<typename T> T Spline<T>::at(float time) const {

    // TODO (Animation): Task 1b
    // Given a time, find the nearest positions & tangent values defined by the control point map.
    if(control_points.size() < 1) return T();//cubic_unit_spline(0.0f, T(), T(), T(), T());
    // Transform them for use with cubic_unit_spline
    if(control_points.size() == 1) return control_points.begin()->second;
    //if query time < start point
    auto end_iter = control_points.end();
    end_iter--;
    auto start_iter = control_points.begin();
    float start_time = start_iter->first;
    float end_time = end_iter->first;

    if(time <= start_time)return start_iter ->second;
    if(time >= end_time)return end_iter->second;
    //when control_points.size() >= 2 and query time between (start,end)
    //k is iterator, point to pair<time, value>
    auto k2 = control_points.upper_bound(time);//upper_bound : >
    if(k2 == control_points.begin())return k2->second;
    if(k2 == control_points.end()) {//k2 == the last element of control_points
        --k2;
        return k2->second;
    }
    auto k1 = k2;//control_points.lower_bound(time);
    k1--;
    //if(k1->first == time)return k1 ->second;

    auto k0 = k1;k0--;
    auto k3 = k2;k3++;

    bool hasK3 = (k2 != end_iter);
    bool hasK0 = (k1 != start_iter);
    float t0 = hasK0 ? ( k0 ->first ):(k1->first - (k2->first - k1->first));
    T p0 = hasK0 ? ( k0 ->second):(k1->second - (k2->second - k1->second));

    //same for k3
    float t3 = hasK3 ? ( k3 -> first):(k2->first + (k2->first - k1->first));
    T p3 = hasK3 ? ( k3 -> second):(k2->second + (k2->second - k1->second));
    ////////General steps
    //Error! t_interval should be normalized
    auto interval = k2->first - k1->first;
    //error! the interval should between k1 and k2, not the whole time interval in animation
    auto m1 = (k2->second - p0)/(k2->first - t0)* interval;//
    auto m2 = (p3 - k1->second)/(t3 - k1->first)* interval;
    //Normalize t
    float t_normalized = (time - k1->first) / interval;//
    return cubic_unit_spline(t_normalized, k1->second, k2->second, m1, m2);
}
