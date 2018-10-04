#include <cmath>
#include "custom_math.h"
using namespace custom_math;

class quaternion
{
public:
	inline quaternion(void) : x(0.0), y(0.0), z(0.0), w(0.0) { /*default constructor*/ }
	inline quaternion(const float &src_x, const float &src_y, const float &src_z, const float &src_w) : x(src_x), y(src_y), z(src_z), w(src_w) { /* custom constructor */ }

	inline float self_dot(void)
	{
		return x*x + y*y + z*z + w*w;
	}

	quaternion operator*(const quaternion &right) const
	{
		quaternion ret;

		ret.x = x*right.x - y*right.y - z*right.z - w*right.w;
		ret.y = x*right.y + y*right.x + z*right.w - w*right.z;
		ret.z = x*right.z - y*right.w + z*right.x + w*right.y;
		ret.w = x*right.w + y*right.z - z*right.y + w*right.x;

		return ret;
	}


	quaternion operator+(const quaternion &right) const
	{
		quaternion ret;

		ret.x = x + right.x;
		ret.y = y + right.y;
		ret.z = z + right.z;
		ret.w = w + right.w;

		return ret;
	}

	float x;
	float y;
	float z;
	float w;
};


quaternion sin(quaternion in)
{
	quaternion ret;

	const float mag_vector = sqrt(in.y*in.y + in.z*in.z + in.w*in.w);

	ret.x = sin(in.x) * cosh(mag_vector);
	ret.y = cos(in.x) * sinh(mag_vector) * in.y / mag_vector;
	ret.z = cos(in.x) * sinh(mag_vector) * in.z / mag_vector;
	ret.w = cos(in.x) * sinh(mag_vector) * in.w / mag_vector;

	return ret;
}


class quaternion_julia_set_params
{
public:
	quaternion_julia_set_params(void)
	{
		z_w = 0;

		c_x = 0.3f;
		c_y = 0.5f;
		c_z = 0.4f;
		c_w = 0.2f;

		x_grid_max = y_grid_max = z_grid_max = 1.5;
		x_grid_min = y_grid_min = z_grid_min = -x_grid_max;
		x_res = y_res = z_res = 150;

		max_iterations = 8;
		threshold = 4;
	}

	float z_w;

	float c_x;
	float c_y;
	float c_z;
	float c_w;

	float x_grid_max, y_grid_max, z_grid_max;
	float x_grid_min, y_grid_min, z_grid_min;
	size_t x_res, y_res, z_res;

	short unsigned int max_iterations;
	float threshold;
};

// Z = Z*Z + C
inline float classic_iter(vector<vector_3> &points, quaternion Z, const quaternion C, const short unsigned int max_iterations, const float threshold)
{
    float length_sq = Z.self_dot(); // In case max_iterations == 0.
    const float threshold_sq = threshold*threshold;
    
    vector_3 p;
    p.x = Z.x;
    p.y = Z.y;
    p.z = Z.z;
    points.push_back(p);
    
    for (short unsigned int i = 0; i < max_iterations; i++)
    {
        Z = Z*Z + C;
    
        vector_3 p;
        p.x = Z.x;
        p.y = Z.y;
        p.z = Z.z;
        points.push_back(p);
        
        // Break (abort) if and when the length of Z squared is larger than threshold squared
        // We deal with squared values in order to save on calls to sqrt
        if ((length_sq = Z.self_dot()) > threshold_sq)
            break;
    }
    
    return sqrt(length_sq);
}

// Z = sin(Z) + C * sin(Z)
inline float sine_plus_iter(vector<vector_3> &points, quaternion Z, const quaternion C, const short unsigned int max_iterations, const float threshold)
{
	float length_sq = Z.self_dot();
	float threshold_sq = threshold*threshold;

    vector_3 p;
    p.x = Z.x;
    p.y = Z.y;
    p.z = Z.z;
    points.push_back(p);
    
	for (short unsigned int i = 0; i < max_iterations; i++)
	{
		Z = sin(Z) + C*sin(Z);
    
        p.x = Z.x;
        p.y = Z.y;
        p.z = Z.z;
        points.push_back(p);
        
		if ((length_sq = Z.self_dot()) > threshold_sq)
			break;
	}

	return sqrt(length_sq);
}


