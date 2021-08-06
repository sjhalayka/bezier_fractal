// Shawn Halayka -- shalayka@gmail.com
// June 26, 2010
//
// This code and data is in the public domain.


#ifndef main_H
#define main_H




#include "bspline_library/Bezier.h"
#include "bspline_library/BSpline.h"
#include "bspline_library/CatmullRom.h"


#include "uv_camera.h"

#include "marching_cubes.h"
using namespace marching_cubes;

#include "eqparse.h"

#include "mesh.h"

#include <fstream>
using std::ofstream;

#include <ios>
using std::ios;

#include <vector>
#include <sstream>
using namespace std;

#include <random>
using std::mt19937;

size_t point_res = 10;

float grid_max = 1.5;

quaternion C(0.2, 0.5, 0.3, 0.4);
unsigned short int max_iterations = 8;
float threshold = 4.0;
float beta = 2.0f;
bool mandelbrot_mode = false;



vector_3 background_colour(1.0, 1.0, 1.0);
float orange[] = { 1.0f, 0.5f, 0.0f, 1.0f };
float mesh_transparent[] = { 0.0f, 0.5f, 1.0f, 0.1f };
float mesh_solid[] = { 0.0f, 0.5f, 1.0f, 1.0f };

float outline_width = 3.0;
static const float outline_colour[] = {0.0, 0.0, 0.0};

bool draw_curves = true;
bool draw_mesh = false;
bool draw_outline = true;
bool draw_axis = true;
bool draw_control_list = true;
bool screenshot_mode = false;

uv_camera main_camera;

GLint win_id = 0;
GLint win_x = 800, win_y = 600;
float camera_w = 4;
float camera_fov = 45;
float camera_x_transform = 0;
float camera_y_transform = 0;
double u_spacer = 0.01;
double v_spacer = 0.5*u_spacer;
double w_spacer = 0.1;
double camera_near = 0.1;
double camera_far = 1000.0;

GLUquadricObj* glu_obj = gluNewQuadric(); // Probably should delete this before app exit... :)

bool lmb_down = false;
bool mmb_down = false;
bool rmb_down = false;
int mouse_x = 0;
int mouse_y = 0;

void idle_func(void);
void init_opengl(const int &width, const int &height);
void reshape_func(int width, int height);
void display_func(void);
void keyboard_func(unsigned char key, int x, int y);
void mouse_func(int button, int state, int x, int y);
void motion_func(int x, int y);
void passive_motion_func(int x, int y);
void draw_objects(bool disable_colouring = false);


vector<triangle> tris;
vector<vertex_3> face_normals;
vector<vertex_3> vertices;
vector<vertex_3> vertex_normals;

vector<vector<vector_4> > all_4d_points;
vector<vector<vector_4> > pos;







float iterate_julia(vector<vector_4>& points, const quaternion& src_Z, const quaternion& src_C, const short unsigned int& max_iterations, const float& threshold, const float exponent)
{
	static quaternion_math qmath;

	quaternion C = src_C;
	quaternion Z = src_Z;

	vector_4 p;
	p.x = Z.x;
	p.y = Z.y;
	p.z = Z.z;
	p.w = Z.w;
	points.push_back(p);

	float len_sq = 0;
	const float threshold_sq = threshold * threshold;

	for (short unsigned int i = 0; i < max_iterations; i++)
	{
		const quaternion P(exponent, 0, 0, 0);

		qmath.pow(&Z, &P, &Z);
		qmath.add(&Z, &C, &Z);

		p.x = Z.x;
		p.y = Z.y;
		p.z = Z.z;
		p.w = Z.w;
		points.push_back(p);

		if ((len_sq = Z.self_dot()) >= threshold_sq)
			break;
	}

	return sqrt(len_sq);
}


float iterate_mandelbrot(vector<vector_4>& points, const quaternion& src_Z, const quaternion &src_C, const short unsigned int& max_iterations, const float& threshold, const float exponent)
{
	static quaternion_math qmath;

	quaternion C = src_Z;
	quaternion Z = quaternion(0, 0, 0, 0);

	vector_4 p;
	p.x = Z.x;
	p.y = Z.y;
	p.z = Z.z;
	p.w = Z.w;
	points.push_back(p);

	float len_sq = 0;
	const float threshold_sq = threshold * threshold;

	for (short unsigned int i = 0; i < max_iterations; i++)
	{
		const quaternion P(exponent, 0, 0, 0);

		qmath.pow(&Z, &P, &Z);
		qmath.add(&Z, &C, &Z);

		p.x = Z.x;
		p.y = Z.y;
		p.z = Z.z;
		p.w = Z.w;
		points.push_back(p);

		if ((len_sq = Z.self_dot()) >= threshold_sq)
			break;
	}

	return sqrt(len_sq);
}

float iterate(bool mandelbrot,
	vector< vector_4 >& trajectory_points,
	quaternion src_Z,
	quaternion src_C,
	const short unsigned int max_iterations,
	const float threshold,
	const float exponent)
{
	trajectory_points.clear();

	if (mandelbrot)
	{
		return iterate_mandelbrot(
			trajectory_points,
			src_Z,
			src_C,
			max_iterations,
			threshold,
			exponent);
	}
	else
	{
		return iterate_julia(
			trajectory_points,
			src_Z,
			src_C,
			max_iterations,
			threshold,
			exponent);
	}
}




// https://stackoverflow.com/questions/785097/how-do-i-implement-a-bézier-curve-in-c
vector_4 getBezierPoint(vector<vector_4> points, float t)
{
	int i = points.size() - 1;

	while (i > 0)
	{
		for (int k = 0; k < i; k++)
		{
			points[k].x += t * (points[k + 1].x - points[k].x);
			points[k].y += t * (points[k + 1].y - points[k].y);
			points[k].z += t * (points[k + 1].z - points[k].z);
			points[k].w += t * (points[k + 1].w - points[k].w);
		}

		i--;
	}

	return points[0];
}



void get_isosurface(
	const bool mandelbrot,
	const float grid_max, 
	const size_t res, 
	const float z_w, 
	const quaternion C, 
	const unsigned short int max_iterations,
	const float threshold,
const float exponent)
{

	const float grid_min = -grid_max;

	const bool make_border = true;

	// When adding a border, use a value that is greater than the threshold.
	const float border_value = 1.0f + threshold;

	tris.clear();
	vector<float> xyplane0(res * res, 0);
	vector<float> xyplane1(res * res, 0);

	const float step_size = (grid_max - grid_min) / (res - 1);

	size_t z = 0;

	quaternion Z(grid_min, grid_min, grid_min, z_w);

	vector<vector_4> points;

	// Calculate xy plane 0.
	for (size_t x = 0; x < res; x++, Z.x += step_size)
	{
		Z.y = grid_min;

		for (size_t y = 0; y < res; y++, Z.y += step_size)
		{
			if (true == make_border && (x == 0 || y == 0 || z == 0 || x == res - 1 || y == res - 1 || z == res - 1))
				xyplane0[x * res + y] = border_value; // 0;
			else
				xyplane0[x * res + y] = iterate(mandelbrot, points, Z, C, max_iterations, threshold, exponent);

			//if (xyplane0[x * res + y] > threshold)
			//	xyplane0[x * res + y] = 0;
			//else
			//	xyplane0[x * res + y] = threshold + 1;

		}
	}

	// Prepare for xy plane 1.
	z++;
	Z.z += step_size;

	size_t box_count = 0;

	// Calculate xy planes 1 and greater.
	for (; z < res; z++, Z.z += step_size)
	{
		Z.x = grid_min;

		cout << "Calculating triangles from xy-plane pair " << z << " of " << res - 1 << endl;

		for (size_t x = 0; x < res; x++, Z.x += step_size)
		{
			Z.y = grid_min;

			for (size_t y = 0; y < res; y++, Z.y += step_size)
			{
				vector<vector_4> points;

				if (true == make_border && (x == 0 || y == 0 || z == 0 || x == res - 1 || y == res - 1 || z == res - 1))
					xyplane1[x * res + y] = border_value; // 0;
				else
					xyplane1[x * res + y] = iterate(mandelbrot, points, Z, C, max_iterations, threshold, exponent);

				//if (xyplane1[x * res + y] > threshold)
				//	xyplane1[x * res + y] = 0;
				//else
				//	xyplane1[x * res + y] = threshold + 1;
			}
		}

		// Calculate triangles for the xy-planes corresponding to z - 1 and z by marching cubes.
		tesselate_adjacent_xy_plane_pair(
			box_count,
			xyplane0, xyplane1,
			z - 1,
			tris,
			threshold, // Use threshold as isovalue.
			grid_min, grid_max, res,
			grid_min, grid_max, res,
			grid_min, grid_max, res);

		// Swap memory pointers (fast) instead of performing a memory copy (slow).
		xyplane1.swap(xyplane0);
	}

	cout << endl;





	if (tris.size() > 0)
	{
		float total_volume = 0;
		float total_area = 0;

		for (size_t i = 0; i < tris.size(); i++)
		{
			vertex_3 a = tris[i].vertex[0];
			vertex_3 b = tris[i].vertex[1];
			vertex_3 c = tris[i].vertex[2];

			total_volume += a.dot(b.cross(c)) / 6.0f;


			vertex_3 a2 = tris[i].vertex[1] - tris[i].vertex[0];
			vertex_3 b2 = tris[i].vertex[2] - tris[i].vertex[0];
			vertex_3 c2 = a2.cross(b2);

			total_area += c2.length() * 0.5f;

		}
		
		cout << "Volume: " << total_volume << endl;
		cout << "Area:   " << total_area << endl;



		get_vertices_and_normals_from_triangles(tris, face_normals, vertices, vertex_normals);
	}

	// Print box-counting dimension
	// Make sure that step_size != 1.0f :)
	//cout << "Box counting dimension: " << logf(static_cast<float>(box_count)) / logf(1.0f / step_size) << endl;
}





void get_points(
	const bool mandelbrot,
	const float grid_max,
	const size_t res,
	const quaternion src_C,
	const unsigned short int max_iterations,
	const float threshold,
	const float exponent)
{
	all_4d_points.clear();
	pos.clear();

	float x_grid_max = grid_max;
	float y_grid_max = grid_max;
	float z_grid_max = grid_max;
	float x_grid_min = -x_grid_max;
	float y_grid_min = -y_grid_max;
	float z_grid_min = -z_grid_max;
	size_t x_res = res;
	size_t y_res = res;
	size_t z_res = res;
	bool make_border = false;

	float z_w = 0;

	get_isosurface(mandelbrot, grid_max, 50, z_w, src_C, max_iterations, threshold, exponent);

	const float x_step_size = (x_grid_max - x_grid_min) / (x_res - 1);
	const float y_step_size = (y_grid_max - y_grid_min) / (y_res - 1);
	const float z_step_size = (z_grid_max - z_grid_min) / (z_res - 1);

	size_t z = 0;

	quaternion Z(x_grid_min, y_grid_min, z_grid_min, z_w);

	for (size_t x = 0; x < x_res; x++, Z.x += x_step_size)
	{
		Z.y = y_grid_min;

		for (size_t y = 0; y < y_res; y++, Z.y += y_step_size)
		{
			vector<vector_4> points;

			float length = iterate(mandelbrot, points, Z, C, max_iterations, threshold, exponent);

			if (length < threshold)
			{
				all_4d_points.push_back(points);
			}
		}
	}

	z++;
	Z.z += z_step_size;

	for (; z < z_res; z++, Z.z += z_step_size)
	{
		Z.x = x_grid_min;

		for (size_t x = 0; x < x_res; x++, Z.x += x_step_size)
		{
			Z.y = y_grid_min;

			for (size_t y = 0; y < y_res; y++, Z.y += y_step_size)
			{
				vector<vector_4> points;

				float length = iterate(mandelbrot, points, Z, C, max_iterations, threshold, exponent);

				if (length < threshold)
				{
					all_4d_points.push_back(points);
				}
			}
		}
	}

	cout << "trajectory count " << all_4d_points.size() << endl;

	size_t orbit_count = 0;

	for (size_t i = 0; i < all_4d_points.size(); i++)
	{
		set<vector_4> point_set;

		for (size_t j = 0; j < all_4d_points[i].size(); j++)
		{
			point_set.insert(all_4d_points[i][j]);
		}

		if (point_set.size() != all_4d_points[i].size())
		{
			cout << point_set.size() << endl;
			orbit_count++;
		}
	}

	cout << "orbit count " << orbit_count << endl;

#define OLD_BEZIER

#ifdef OLD_BEZIER



	for (size_t i = 0; i < all_4d_points.size(); i++)
	{
		vector<vector_4> p;

		for (float t = 0; t <= 1.0f; t += 0.01f)
		//for (float t = 0; t <= 0.85f; t += 0.01f)
		{
			vector_4 v = getBezierPoint(all_4d_points[i], t);
			p.push_back(v);
		}

		pos.push_back(p);
	}


#else


	for (size_t i = 0; i < all_4d_points.size(); i++)
	{
		Curve* curve = new CatmullRom();
		curve->set_steps(10);

		for (size_t j = 0; j < all_4d_points[i].size(); j++)
			curve->add_way_point(Vector(all_4d_points[i][j].x, all_4d_points[i][j].y, all_4d_points[i][j].z));

//		if (is_cycle[i])
///		{
//			for (size_t j = 0; j < all_4d_points[i].size(); j++)
//				curve->add_way_point(Vector(all_4d_points[i][j].x, all_4d_points[i][j].y, all_4d_points[i][j].z));
//		}

		vector<vector_4> p;

		for (int i = 0; i < curve->node_count(); i++)
			p.push_back(vector_4(curve->node(i).x, curve->node(i).y, curve->node(i).z, 0));

		if (p.size() == 0)
			p.push_back(all_4d_points[i][0]);

		pos.push_back(p);

		delete curve;
	}

#endif
}


// TODO: fix camera bug where portrait mode crashes.
void take_screenshot(size_t num_cams_wide, const char *filename, const bool reverse_rows = false)
{
	screenshot_mode = true;

	

	// Set up Targa TGA image data.
	unsigned char  idlength = 0;
	unsigned char  colourmaptype = 0;
	unsigned char  datatypecode = 2;
	unsigned short int colourmaporigin = 0;
	unsigned short int colourmaplength = 0;
	unsigned char  colourmapdepth = 0;
	unsigned short int x_origin = 0;
	unsigned short int y_origin = 0;

	cout << "Image size: " << static_cast<size_t>(win_x)*num_cams_wide << "x" << static_cast<size_t>(win_y)*num_cams_wide << " pixels" << endl;

	if (static_cast<size_t>(win_x)*num_cams_wide > static_cast<unsigned short>(-1) ||
		static_cast<size_t>(win_y)*num_cams_wide > static_cast<unsigned short>(-1))
	{
		cout << "Image too large. Maximum width and height is " << static_cast<unsigned short>(-1) << endl;
		return;
	}

	unsigned short int px = win_x*static_cast<unsigned short>(num_cams_wide);
	unsigned short int py = win_y*static_cast<unsigned short>(num_cams_wide);
	unsigned char  bitsperpixel = 24;
	unsigned char  imagedescriptor = 0;
	vector<char> idstring;

	size_t num_bytes = 3*px*py;
	vector<unsigned char> pixel_data(num_bytes);

	// Adjust some parameters for large screen format.
	bool temp_draw_control_list = draw_control_list;
	draw_control_list = false;

	float temp_outline_width = outline_width;
	outline_width = 6;

	vector<unsigned char> fbpixels(3*win_x*win_y);

	const size_t total_cams = num_cams_wide * num_cams_wide;
	size_t cam_count = 0;
	// Loop through subcameras.
	for(size_t cam_num_x = 0; cam_num_x < num_cams_wide; cam_num_x++)
	{
		for(size_t cam_num_y = 0; cam_num_y < num_cams_wide; cam_num_y++)
		{
			cout << "Camera: " << cam_count + 1 << " of " << total_cams << endl;

			// Set up camera, draw, then copy the frame buffer.
			main_camera.Set_Large_Screenshot(num_cams_wide, cam_num_x, cam_num_y);
			display_func();
			glReadPixels(0, 0, win_x, win_y, GL_RGB, GL_UNSIGNED_BYTE, &fbpixels[0]);

			// Copy pixels to large image.
			for(GLint i = 0; i < win_x; i++)
			{
				for(GLint j = 0; j < win_y; j++)
				{
					size_t fb_index = 3*(j*win_x + i);

					size_t screenshot_x = cam_num_x*win_x + i;
					size_t screenshot_y = cam_num_y*win_y + j;
					size_t screenshot_index = 3*(screenshot_y*(win_x*num_cams_wide) + screenshot_x);

					pixel_data[screenshot_index] = fbpixels[fb_index + 2];
					pixel_data[screenshot_index + 1] = fbpixels[fb_index + 1];
					pixel_data[screenshot_index + 2] = fbpixels[fb_index ];
				}
			}

			cam_count++;
		}

	}

	screenshot_mode = false;

	// Restore the parameters.
	draw_control_list = temp_draw_control_list;
	outline_width = temp_outline_width;
	main_camera.Set();

	// Write Targa TGA file to disk.
	ofstream out(filename, ios::binary);

	if(!out.is_open())
	{
		cout << "Failed to open TGA file for writing: " << filename << endl;
		return;
	}

	out.write(reinterpret_cast<char *>(&idlength), 1);
	out.write(reinterpret_cast<char *>(&colourmaptype), 1);
	out.write(reinterpret_cast<char *>(&datatypecode), 1);
	out.write(reinterpret_cast<char *>(&colourmaporigin), 2);
	out.write(reinterpret_cast<char *>(&colourmaplength), 2);
	out.write(reinterpret_cast<char *>(&colourmapdepth), 1);
	out.write(reinterpret_cast<char *>(&x_origin), 2);
	out.write(reinterpret_cast<char *>(&y_origin), 2);
	out.write(reinterpret_cast<char *>(&px), 2);
	out.write(reinterpret_cast<char *>(&py), 2);
	out.write(reinterpret_cast<char *>(&bitsperpixel), 1);
	out.write(reinterpret_cast<char *>(&imagedescriptor), 1);

	out.write(reinterpret_cast<char *>(&pixel_data[0]), num_bytes);


}

void idle_func(void)
{
	glutPostRedisplay();
}




void init_opengl(const int &width, const int &height)
{
	win_x = width;
	win_y = height;

	if(win_x < 1)
		win_x = 1;

	if(win_y < 1)
		win_y = 1;

	glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_ALPHA|GLUT_DEPTH);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(win_x, win_y);
	win_id = glutCreateWindow("2D v4");

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glDepthMask(GL_TRUE);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_POLYGON_SMOOTH);

	float light_colour[] = {0.9f, 0.9f, 0.9f, 1.0f};
	float light0_position[] = {0.0f, 0.0f, 1.0f, 0.0f };
	float light1_position[] = {0.0f, 0.0f, -1.0f, 0.0f };
	float light2_position[] = {1.0f, 0.0f, 0.0f, 0.0f };
	float light3_position[] = {-1.0f, 0.0f, 0.0f, 0.0f };

	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_colour);
	glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_colour);
	glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, light_colour);
	glLightfv(GL_LIGHT3, GL_POSITION, light3_position);
	glLightfv(GL_LIGHT3, GL_DIFFUSE, light_colour);

	float light_colour2[] = {0.5f, 0.5f, 0.5f, 1.0f};
	float light4_position[] = {0.0f, 1.0f, 0.0f, 0.0f };
	float light5_position[] = {0.0f, -1.0f, 0.0f, 0.0f };

	glLightfv(GL_LIGHT4, GL_POSITION, light4_position);
	glLightfv(GL_LIGHT4, GL_DIFFUSE, light_colour2);
	glLightfv(GL_LIGHT5, GL_POSITION, light5_position);
	glLightfv(GL_LIGHT5, GL_DIFFUSE, light_colour2);

	glClearColor(static_cast<float>(background_colour.x), static_cast<float>(background_colour.y), static_cast<float>(background_colour.z), 1.0f);
	glClearDepth(1.0f);

	main_camera.Set(0, 0, camera_w, camera_fov, win_x, win_y, camera_near, camera_far);

	get_points(mandelbrot_mode, grid_max, point_res, C, max_iterations, threshold, beta);
}

void reshape_func(int width, int height)
{
	win_x = width;
	win_y = height;

	if(win_x < 1)
		win_x = 1;

	if(win_y < 1)
		win_y = 1;

	glutSetWindow(win_id);
	glutReshapeWindow(win_x, win_y);
	glViewport(0, 0, win_x, win_y);

	main_camera.Set(main_camera.u, main_camera.v, main_camera.w, main_camera.fov, win_x, win_y, camera_near, camera_far);
}

// Text drawing code originally from "GLUT Tutorial -- Bitmap Fonts and Orthogonal Projections" by A R Fernandes
void render_string(int x, const int y, void *font, const string &text)
{
	for(size_t i = 0; i < text.length(); i++)
	{
		glRasterPos2i(x, y);
		glutBitmapCharacter(font, text[i]);
		x += glutBitmapWidth(font, text[i]) + 1;
	}
}
// End text drawing code.

void display_func(void)
{
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);


	if(true == draw_outline)
	{
		glDisable(GL_DEPTH_TEST);

		// Draw outline code from NeHe lesson 37:
		// http://nehe.gamedev.net/data/lessons/lesson.asp?lesson=37
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glLineWidth(outline_width);
		glCullFace(GL_BACK);
		glPolygonMode(GL_FRONT, GL_LINE);
		glColor3fv(&outline_colour[0]);

		draw_objects(true);


		glPopAttrib();
		// End draw outline code.

		glEnable(GL_DEPTH_TEST);
	}


	// Draw the model's components using OpenGL/GLUT primitives.



	draw_objects();




	if(true == draw_outline)
	{
		// Draw outline code from NeHe lesson 37:
		// http://nehe.gamedev.net/data/lessons/lesson.asp?lesson=37
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glLineWidth(outline_width);
		glCullFace(GL_FRONT);
		glPolygonMode(GL_BACK, GL_LINE);

		draw_objects(true);

		glPopAttrib();
		// End draw outline code.
	}

	if(true == draw_control_list)
	{
		// Text drawing code originally from "GLUT Tutorial -- Bitmap Fonts and Orthogonal Projections" by A R Fernandes
		// http://www.lighthouse3d.com/opengl/glut/index.php?bmpfontortho
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluOrtho2D(0, win_x, 0, win_y);
		glScalef(1, -1, 1); // Neat. :)
		glTranslatef(0, -static_cast<float>(win_y), 0); // Neat. :)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		glColor3f(0, 0, 0);

		size_t break_size = 13;
		size_t start = 20;

		render_string(10, start, GLUT_BITMAP_HELVETICA_12, string("Keyboard controls:"));

		render_string(10, start + 1*break_size, GLUT_BITMAP_HELVETICA_10, string("H: Draw outlines"));
		render_string(10, start + 2*break_size, GLUT_BITMAP_HELVETICA_10, string("J: Draw axis"));
		render_string(10, start + 3*break_size, GLUT_BITMAP_HELVETICA_10, string("K: Draw this list"));
		render_string(10, start + 4*break_size, GLUT_BITMAP_HELVETICA_10, string("L: Take screenshot"));





		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		//glMatrixMode(GL_MODELVIEW);
		// End text drawing code.
	}


	glEnable(GL_DEPTH_TEST);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHT2);
	glEnable(GL_LIGHT3);
	glEnable(GL_LIGHT4);
	glEnable(GL_LIGHT5);

	glEnable(GL_ALPHA);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	if(draw_mesh)
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mesh_solid);
	else
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mesh_transparent);

	glBegin(GL_TRIANGLES);

	for (size_t i = 0; i < tris.size(); i++)
	{
		size_t v_index0 = tris[i].vertex[0].index;
		size_t v_index1 = tris[i].vertex[1].index;
		size_t v_index2 = tris[i].vertex[2].index;

		glNormal3f(vertex_normals[v_index0].x, vertex_normals[v_index0].y, vertex_normals[v_index0].z);
		glVertex3f(vertices[v_index0].x, vertices[v_index0].y, vertices[v_index0].z);
		glNormal3f(vertex_normals[v_index1].x, vertex_normals[v_index1].y, vertex_normals[v_index1].z);
		glVertex3f(vertices[v_index1].x, vertices[v_index1].y, vertices[v_index1].z);
		glNormal3f(vertex_normals[v_index2].x, vertex_normals[v_index2].y, vertex_normals[v_index2].z);
		glVertex3f(vertices[v_index2].x, vertices[v_index2].y, vertices[v_index2].z);
	}

	glEnd();





	if (false == screenshot_mode)
	{
		glFlush();
		glutSwapBuffers();
	}
}

void keyboard_func(unsigned char key, int x, int y)
{
	switch(tolower(key))
	{
	case 'y':
	{
		draw_mesh = !draw_mesh;
		break;
	}
	case 'u':
	{
		draw_curves = !draw_curves;
		break;
	}
	case 'h':
		{
			draw_outline = !draw_outline;
			break;
		}
	case 'j':
		{
			draw_axis = !draw_axis;
			break;
		}
	case 'k':
		{
			draw_control_list = !draw_control_list;
			break;
		}
	case 'n':
	{
		get_points(mandelbrot_mode, grid_max, point_res, C, max_iterations, threshold, beta);
		take_screenshot(8, "screenshot.tga");
		get_points(mandelbrot_mode, grid_max, point_res, C, max_iterations, threshold, beta);
		break;
	}
	case 'm':
		{
			get_points(mandelbrot_mode, grid_max, 20, C, max_iterations, threshold, beta);
			take_screenshot(8, "screenshot.tga");
			get_points(mandelbrot_mode, grid_max, point_res, C, max_iterations, threshold, beta);
			break;
		}
	case 'b':
	{
		vector<float> betas = { -8, -7, -6, -5, -4, -3, 2, 3, 4, 5, 6, 7, 8 };

		for (size_t i = 0; i < betas.size(); i++)
		{
			ostringstream oss;
			oss << betas[i];

			string filename = "julia3d_" + oss.str() + ".tga";

			cout << filename << endl;

			get_points(false, grid_max, 20, C, max_iterations, threshold, betas[i]);
			take_screenshot(2, filename.c_str());
		}

		for (size_t i = 0; i < betas.size(); i++)
		{
			ostringstream oss;
			oss << betas[i];

			string filename = "mandelbrot3d_" + oss.str() + ".tga";

			cout << filename << endl;

			get_points(true, grid_max, 20, C, max_iterations, threshold, betas[i]);
			take_screenshot(2, filename.c_str());
		}

		break;
	}
	default:
		break;
	}
}

void mouse_func(int button, int state, int x, int y)
{
	if(GLUT_LEFT_BUTTON == button)
	{
		if(GLUT_DOWN == state)
			lmb_down = true;
		else
			lmb_down = false;
	}
	else if(GLUT_MIDDLE_BUTTON == button)
	{
		if(GLUT_DOWN == state)
			mmb_down = true;
		else
			mmb_down = false;
	}
	else if(GLUT_RIGHT_BUTTON == button)
	{
		if(GLUT_DOWN == state)
			rmb_down = true;
		else
			rmb_down = false;
	}
}

void motion_func(int x, int y)
{
	int prev_mouse_x = mouse_x;
	int prev_mouse_y = mouse_y;

	mouse_x = x;
	mouse_y = y;

	int mouse_delta_x = mouse_x - prev_mouse_x;
	int mouse_delta_y = prev_mouse_y - mouse_y;

	if(true == lmb_down && (0 != mouse_delta_x || 0 != mouse_delta_y))
	{
		main_camera.u -= static_cast<float>(mouse_delta_y)*u_spacer;
		main_camera.v += static_cast<float>(mouse_delta_x)*v_spacer;
	}
	else if(true == rmb_down && (0 != mouse_delta_y))
	{
		main_camera.w -= static_cast<float>(mouse_delta_y)*w_spacer;
	}

	main_camera.Set(); // Calculate new camera vectors.
}

void passive_motion_func(int x, int y)
{
	mouse_x = x;
	mouse_y = y;
}



class RGB
{
public:
	unsigned char r, g, b;
};

RGB HSBtoRGB(unsigned short int hue_degree, unsigned char sat_percent, unsigned char bri_percent)
{
	float R = 0.0f;
	float G = 0.0f;
	float B = 0.0f;

	if (hue_degree > 359)
		hue_degree = 359;

	if (sat_percent > 100)
		sat_percent = 100;

	if (bri_percent > 100)
		bri_percent = 100;

	float hue_pos = 6.0f - ((static_cast<float>(hue_degree) / 359.0f) * 6.0f);

	if (hue_pos >= 0.0f && hue_pos < 1.0f)
	{
		R = 255.0f;
		G = 0.0f;
		B = 255.0f * hue_pos;
	}
	else if (hue_pos >= 1.0f && hue_pos < 2.0f)
	{
		hue_pos -= 1.0f;

		R = 255.0f - (255.0f * hue_pos);
		G = 0.0f;
		B = 255.0f;
	}
	else if (hue_pos >= 2.0f && hue_pos < 3.0f)
	{
		hue_pos -= 2.0f;

		R = 0.0f;
		G = 255.0f * hue_pos;
		B = 255.0f;
	}
	else if (hue_pos >= 3.0f && hue_pos < 4.0f)
	{
		hue_pos -= 3.0f;

		R = 0.0f;
		G = 255.0f;
		B = 255.0f - (255.0f * hue_pos);
	}
	else if (hue_pos >= 4.0f && hue_pos < 5.0f)
	{
		hue_pos -= 4.0f;

		R = 255.0f * hue_pos;
		G = 255.0f;
		B = 0.0f;
	}
	else
	{
		hue_pos -= 5.0f;

		R = 255.0f;
		G = 255.0f - (255.0f * hue_pos);
		B = 0.0f;
	}

	if (100 != sat_percent)
	{
		if (0 == sat_percent)
		{
			R = 255.0f;
			G = 255.0f;
			B = 255.0f;
		}
		else
		{
			if (255.0f != R)
				R += ((255.0f - R) / 100.0f) * (100.0f - sat_percent);
			if (255.0f != G)
				G += ((255.0f - G) / 100.0f) * (100.0f - sat_percent);
			if (255.0f != B)
				B += ((255.0f - B) / 100.0f) * (100.0f - sat_percent);
		}
	}

	if (100 != bri_percent)
	{
		if (0 == bri_percent)
		{
			R = 0.0f;
			G = 0.0f;
			B = 0.0f;
		}
		else
		{
			if (0.0f != R)
				R *= static_cast<float>(bri_percent) / 100.0f;
			if (0.0f != G)
				G *= static_cast<float>(bri_percent) / 100.0f;
			if (0.0f != B)
				B *= static_cast<float>(bri_percent) / 100.0f;
		}
	}

	if (R < 0.0f)
		R = 0.0f;
	else if (R > 255.0f)
		R = 255.0f;

	if (G < 0.0f)
		G = 0.0f;
	else if (G > 255.0f)
		G = 255.0f;

	if (B < 0.0f)
		B = 0.0f;
	else if (B > 255.0f)
		B = 255.0f;

	RGB rgb;

	rgb.r = static_cast<unsigned char>(R);
	rgb.g = static_cast<unsigned char>(G);
	rgb.b = static_cast<unsigned char>(B);

	return rgb;
}




// This render mode won't apply to a curved 3D space.
void draw_objects(bool disable_colouring)
{





	if (false == disable_colouring)
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);
		glEnable(GL_LIGHT2);
		glEnable(GL_LIGHT3);
		glEnable(GL_LIGHT4);
		glEnable(GL_LIGHT5);
	}
	else
	{
		glColor3f(0.0f, 0.0f, 0.0f);
		glDisable(GL_LIGHTING);
	}

	static const float rad_to_deg = 180.0f / static_cast<float>(pi);

	glPushMatrix();

	glTranslatef(camera_x_transform, camera_y_transform, 0);



	//for (size_t i = 0; i < grid_vertices.size(); i++)
	//{
	//	glPushMatrix();
	//	glTranslatef(grid_vertices[i].x, grid_vertices[i].y, grid_vertices[i].z);
	//	glutSolidSphere(0.025, 16, 16);
	//	glPopMatrix();
	//}



	if (draw_curves)
	{
		if (false == disable_colouring)
			glMaterialfv(GL_FRONT, GL_DIFFUSE, orange);

		for (size_t i = 0; i < pos.size(); i++)
		{
			for (size_t j = 0; j < pos[i].size() - 1; j++)
			{
				double t = j / static_cast<double>(pos[i].size() - 1);

				//set<vector_4> point_set;

				//for (size_t j = 0; j < all_4d_points[i].size(); j++)
				//{
				//	point_set.insert(all_4d_points[i][j]);
				//}

				//if (point_set.size() == all_4d_points[i].size())
				//	continue;

				// double t = static_cast<float>(point_set.size()) / static_cast<float>(all_4d_points[i].size());

				RGB rgb = HSBtoRGB(static_cast<unsigned short>(300.f * t), 75, 100);

				float colour[] = { rgb.r / 255.0f, rgb.g / 255.0f, rgb.b / 255.0f, 1.0f };

				glMaterialfv(GL_FRONT, GL_DIFFUSE, colour);

				vector_4 line = pos[i][j + 1] - pos[i][j];

				glPushMatrix();
				glTranslatef(static_cast<float>(pos[i][j].x), static_cast<float>(pos[i][j].y), static_cast<float>(pos[i][j].z));

				float line_len = static_cast<float>(line.length());
				line.normalize();

				float yaw = 0.0f;

				if (fabsf(static_cast<float>(line.x)) < 0.00001f && fabsf(static_cast<float>(line.z)) < 0.00001f)
					yaw = 0.0f;
				else
					yaw = atan2f(static_cast<float>(line.x), static_cast<float>(line.z));

				float pitch = -atan2f(static_cast<float>(line.y), static_cast<float>(sqrt(line.x * line.x + line.z * line.z)));

				glRotatef(yaw * rad_to_deg, 0.0f, 1.0f, 0.0f);
				glRotatef(pitch * rad_to_deg, 1.0f, 0.0f, 0.0f);

				if (j == 0)
					glutSolidSphere(0.005 * 1.5, 16, 16);

				if (j < pos[i].size() - 2)
					gluCylinder(glu_obj, 0.005, 0.005, line_len, 20, 2);
				else
					glutSolidCone(0.005 * 4, 0.005 * 8, 20, 20);

				glPopMatrix();
			}

		}


		//for (size_t i = 0; i < all_4d_points.size(); i++)
		//{
		//	for (size_t j = 0; j < all_4d_points[i].size() - 1; j++)
		//	{
		//		double t = j / static_cast<double>(all_4d_points[i].size() - 1);

		//		RGB rgb = HSBtoRGB(static_cast<unsigned short>(300.f * t), 75, 100);

		//		float colour[] = { rgb.r / 255.0f, rgb.g / 255.0f, rgb.b / 255.0f, 1.0f };

		//		glMaterialfv(GL_FRONT, GL_DIFFUSE, colour);

		//		vector_4 line = all_4d_points[i][j + 1] - all_4d_points[i][j];

		//		glPushMatrix();
		//		glTranslatef(static_cast<float>(all_4d_points[i][j].x), static_cast<float>(all_4d_points[i][j].y), static_cast<float>(all_4d_points[i][j].z));

		//		float line_len = static_cast<float>(line.length());
		//		line.normalize();

		//		float yaw = 0.0f;

		//		if (fabsf(static_cast<float>(line.x)) < 0.00001f && fabsf(static_cast<float>(line.z)) < 0.00001f)
		//			yaw = 0.0f;
		//		else
		//			yaw = atan2f(static_cast<float>(line.x), static_cast<float>(line.z));

		//		float pitch = -atan2f(static_cast<float>(line.y), static_cast<float>(sqrt(line.x * line.x + line.z * line.z)));

		//		glRotatef(yaw * rad_to_deg, 0.0f, 1.0f, 0.0f);
		//		glRotatef(pitch * rad_to_deg, 1.0f, 0.0f, 0.0f);

		//		if (j == 0)
		//			glutSolidSphere(0.005 * 1.5, 16, 16);

		//		if (j < all_4d_points[i].size() - 2)
		//			gluCylinder(glu_obj, 0.005, 0.005, line_len, 20, 2);
		//		else
		//			glutSolidCone(0.005 * 4, 0.005 * 8, 20, 20);

		//		glPopMatrix();
		//	}

		//}

	}


	if (draw_mesh)
	{
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mesh_solid);

		glBegin(GL_TRIANGLES);

		for (size_t i = 0; i < tris.size(); i++)
		{
			size_t v_index0 = tris[i].vertex[0].index;
			size_t v_index1 = tris[i].vertex[1].index;
			size_t v_index2 = tris[i].vertex[2].index;

			glNormal3f(vertex_normals[v_index0].x, vertex_normals[v_index0].y, vertex_normals[v_index0].z);
			glVertex3f(vertices[v_index0].x, vertices[v_index0].y, vertices[v_index0].z);
			glNormal3f(vertex_normals[v_index1].x, vertex_normals[v_index1].y, vertex_normals[v_index1].z);
			glVertex3f(vertices[v_index1].x, vertices[v_index1].y, vertices[v_index1].z);
			glNormal3f(vertex_normals[v_index2].x, vertex_normals[v_index2].y, vertex_normals[v_index2].z);
			glVertex3f(vertices[v_index2].x, vertices[v_index2].y, vertices[v_index2].z);
		}

		glEnd();
	}







	glDisable(GL_LIGHTING);

	// If we do draw the axis at all, make sure not to draw its outline.
	if(draw_axis && false == disable_colouring)
	{
		glEnable(GL_ALPHA);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glLineWidth(outline_width);

		glBegin(GL_LINES);

		//glColor4f(0, 0, 0, 0.5);

		glColor4f(1, 0, 0, 0.25);
		glVertex3f(0, 0, 0);
		glVertex3f(1, 0, 0);
		glColor4f(0, 1, 0, 0.25);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 1, 0);
		glColor4f(0, 0, 1, 0.25);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, 1);

		//glColor4f(0, 0, 0, 0.25);
		//glVertex3f(0, 0, 0);
		//glVertex3f(-1, 0, 0);
		//glVertex3f(0, 0, 0);
		//glVertex3f(0, -1, 0);
		//glVertex3f(0, 0, 0);
		//glVertex3f(0, 0, -1);

		glEnd();

		glDisable(GL_BLEND);
		glDisable(GL_ALPHA);
	}

	glPopMatrix();
}








#endif