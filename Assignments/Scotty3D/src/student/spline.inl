
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

    // Given a time, find the nearest positions & tangent values
    // defined by the control point map.
    if(control_points.size() == 0) return T();//cubic_unit_spline(0.0f, T(), T(), T(), T());
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
    time = clamp(time, start_time, end_time);//test
    //when control_points.size() >= 2 and query time between (start,end)
    //k is iterator, point to pair<time, value>
    auto k2 = control_points.upper_bound(time);
    if(k2 == control_points.end()) --k2;//test
    auto k1 = k2;//control_points.lower_bound(time);
    k1--;
    if(k1->first == time)return k1 ->second;

    auto k0 = k1;
    auto k3 = k2;

    bool hasK3 = (k3 != end_iter);//(control_points.upper_bound(k2->first) != control_points.end());
    bool hasK0 = (k0 != start_iter);//(control_points.lower_bound(k1->first) != control_points.end());

    //if k0 exists: k0 = k1-1

    float t0 = hasK0 ? ( (--k0) ->first ):(k1->first - (k2->first - k1->first));
    T p0 = hasK0 ? ( k0 ->second):(k1->second - (k2->second - k1->second));
    // if(control_points.lower_bound(k1->first) == control_points.end()){
    //     t0 = k1->first - (k2->first - k1->first);
    //     p0 = k1->second - (k2->second - k1->second);
    // }
    // //has k0
    // else{
    //     auto k0 = control_points.lower_bound(k1->first);
    //     t0 = k0->first;
    //     p0 = k0->second;
    // }
    //same for k3

    float t3 = hasK3 ? ( (++k3) -> first):(k2->first + (k2->first - k1->first));
    T p3 = hasK3 ? ( k3 -> second):(k2->second + (k2->second - k1->second));
    // if(control_points.upper_bound(k2->first) == control_points.end()){
    //     t3 = k2->first + (k2->first - k1->first);
    //     p3 = k2->second + (k2->second - k1->second);
    // }
    // else{
    //     auto k3 = control_points.upper_bound(k2->first);
    //     t3 = k3->first;
    //     p3 = k3->second;
    // }
    ////////General steps
    //Error! t_interval should be normalized
    auto m1 = (k2->second - p0)/(k2->first - t0)* (end_time - start_time);//* (end_time - start_time)
    auto m2 = (p3 - k1->second)/(t3 - k1->first)* (end_time - start_time);
    //Normalize t
    float t_normalized = (time - start_time) / (end_time - start_time);// Error! time / (end_time - start_time)
    return cubic_unit_spline(t_normalized, k1->second, k2->second, m1, m2);
}
