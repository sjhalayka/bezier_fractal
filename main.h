// Shawn Halayka -- shalayka@gmail.com
// June 26, 2010
//
// This code and data is in the public domain.


#ifndef main_H
#define main_H


#include "uv_camera.h"
//#include "qjs.h"
#include "eqparse.h"

#include <fstream>
using std::ofstream;

#include <ios>
using std::ios;

#include <vector>
#include <sstream>
using namespace std;

#include <opencv2/opencv.hpp>
using namespace cv;
#pragma comment(lib, "opencv_world340.lib")


vector_3 background_colour(1.0, 1.0, 1.0);
float sphere_colour[] = { 1.0f, 0.5f, 0.0f, 1.0f };

float outline_width = 3.0;
static const float outline_colour[] = {0.0, 0.0, 0.0};

bool draw_outline = true;
bool draw_axis = true;
bool draw_control_list = true;


uv_camera main_camera;

GLint win_id = 0;
GLint win_x = 800, win_y = 600;
float camera_w = 3;
float camera_fov = 45;
float camera_x_transform = 0;
float camera_y_transform = 0;
double u_spacer = 0.01;
double v_spacer = 0.5*u_spacer;
double w_spacer = 0.1;
double camera_near = 0.1;
double camera_far = 10.0;

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



vector<vector<vector_4> > all_4d_points;
vector<vector<vector_4> > pos;

vector<float> total_lengths;
vector<float> total_distances;
vector<float> dist_by_len;
vector<float> mags;

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



float mean(const vector<float> &src)
{
	float m = 0;
	float size = static_cast<float>(src.size());

	for (size_t i = 0; i < src.size(); i++)
		m += src[i];

	m /= size;

	return m;
}

float standard_deviation(const vector<float> &src)
{
	float m = mean(src);
	float size = static_cast<float>(src.size());

	float sq_diff = 0;

	for (size_t i = 0; i < src.size(); i++)
	{
		float diff = src[i] - m;
		sq_diff += diff*diff;
	}

	sq_diff /= size;

	return sqrtf(sq_diff);
}

float variance(const vector<float> &src)
{
	float s = standard_deviation(src);

	return s*s;
}

// https://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
float skewness(const vector<float> &src)
{
	float m = mean(src);
	float s = standard_deviation(src);
	float size = static_cast<float>(src.size());

	float cube_diff = 0;

	for (size_t i = 0; i < src.size(); i++)
	{
		float diff = src[i] - m;
		cube_diff += diff*diff*diff;
	}

	cube_diff /= size;
	cube_diff /= s*s*s;

	return cube_diff;
}

// https://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
float kurtosis(const vector<float> &src)
{
	float m = mean(src);
	float s = standard_deviation(src);
	float size = static_cast<float>(src.size());

	float fourth_diff = 0;

	for (size_t i = 0; i < src.size(); i++)
	{
		float diff = src[i] - m;
		fourth_diff += diff*diff*diff*diff;
	}

	fourth_diff /= size;
	fourth_diff /= s*s*s*s;

	fourth_diff -= 3.0f;

	return fourth_diff;
}



void get_points(size_t res)
{
	all_4d_points.clear();
	pos.clear();

	total_lengths.clear();
	total_distances.clear();
	dist_by_len.clear();

	float x_grid_max = 1.5;
	float y_grid_max = 1.5;
	float z_grid_max = 1.5;
	float x_grid_min = -x_grid_max;
	float y_grid_min = -y_grid_max;
	float z_grid_min = -z_grid_max;
	size_t x_res = res;
	size_t y_res = res;
	size_t z_res = res;
	bool make_border = true;

	float z_w = 0;
	quaternion C;
	C.x = 0.3f;
	C.y = 0.5f;
	C.z = 0.4f;
	C.w = 0.2f;
	unsigned short int max_iterations = 8;
	float threshold = 4;


	string error_string;
	quaternion_julia_set_equation_parser eqparser;
	if (false == eqparser.setup("Z = Z^5 + C", error_string, C))
	{
		cout << "Equation error: " << error_string << endl;
		return;
	}

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

			float length = eqparser.iterate(points, Z, max_iterations, threshold);

			if (length < threshold)
			{
				mags.push_back(length);
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

				float length = eqparser.iterate(points, Z, max_iterations, threshold);

				if (length < threshold)
				{
					mags.push_back(length);
					all_4d_points.push_back(points);
				}
			}
		}
	}




	//for (size_t i = 0; i < all_4d_points.size(); i++)
	//{
	//	vector<vector_4> p;

	//	for (float t = 0; t <= 1.0; t += 0.01)
	//	{
	//		vector_4 v = getBezierPoint(all_4d_points[i], t);
	//		p.push_back(v);
	//	}

	//	pos.push_back(p);
	//}

	pos = all_4d_points;


	//for (size_t i = 0; i < pos.size(); i++)
	//{
	//	for (size_t j = 0; j < pos[i].size(); j++)
	//	{
	//		vector_4 mapped;
	//		float x = pos[i][j].x;
	//		float y = pos[i][j].y;
	//		float z = pos[i][j].z;
	//		float w = pos[i][j].w;

	//		mapped.x = 2.0 * (x*y + z*w);
	//		mapped.y = 2.0 * (x*w + y*z);
	//		mapped.z = (x*x + z*z) - (y*y + w*w);
	//		mapped.w = 0;

	//		pos[i][j] = mapped;
	//	}

	//}

	for (size_t i = 0; i < pos.size(); i++)
	{
		float total_len = 0;

		for (size_t j = 0; j < pos[i].size() - 1; j++)
		{
			vector_4 line = pos[i][j + 1] - pos[i][j];

			float line_len = line.length();		
			total_len += line_len;
		}

		total_lengths.push_back(total_len);

		float distance = (pos[pos.size() - 1][0] - pos[i][0]).length();
		total_distances.push_back(distance);
	}

	for (size_t i = 0; i < total_lengths.size(); i++)
		dist_by_len.push_back(total_distances[i] / total_lengths[i]);

	total_lengths = mags;
//	total_lengths = total_distances;
//	total_lengths = dist_by_len;

	cout << endl;
	cout << "mean:     " << mean(total_lengths) << endl;
	cout << "std dev:  " << standard_deviation(total_lengths) << endl;
	cout << "variance: " << variance(total_lengths) << endl;
	cout << "skewness: " << skewness(total_lengths) << endl;
	cout << "kurtosis: " << kurtosis(total_lengths) << endl;







	float max_length = 0;

	for (size_t i = 0; i < total_lengths.size(); i++)
	{
		if (total_lengths[i] > max_length)
			max_length = total_lengths[i];
	}

	for (size_t i = 0; i < total_lengths.size(); i++)
	{
		total_lengths[i] /= max_length;
		total_lengths[i] *= 255.0f;
		total_lengths[i] = floorf(total_lengths[i]);
	}





	Mat data(1, total_lengths.size(), CV_8UC1, Scalar(0));

	for (size_t i = 0; i < total_lengths.size(); i++)
		data.at<unsigned char>(0, i) = static_cast<unsigned char>(total_lengths[i]);



	int histSize = 256;
	float range[] = { 0, 256 }; //the upper boundary is exclusive
	const float* histRange = { range };
	bool uniform = true, accumulate = false;
	Mat hist;
	calcHist(&data, 1, 0, Mat(), hist, 1, &histSize, &histRange, uniform, accumulate);

	int hist_w = 600, hist_h = 600;
	int bin_w = cvRound((double)hist_w / histSize);
	Mat histImage(hist_h, hist_w, CV_8UC3, Scalar(255, 255, 255));
	normalize(hist, hist, 0, histImage.rows, NORM_MINMAX, -1, Mat());

	float largest_hist = 0;
	float largest_hist_j = 0;

	for (int j = 0; j < hist.rows; j++)
	{
		for (int i = 0; i < hist.cols; i++)
		{
			if (hist.at<float>(j, i) > largest_hist)
			{
				largest_hist = hist.at<float>(j, i);
				largest_hist_j = j;
			}
		}
	}

	for (int i = 1; i < histSize; i++)
	{
		line(histImage, Point(bin_w*(i - 1), hist_h - cvRound(hist.at<float>(i - 1))),
			Point(bin_w*(i), hist_h - cvRound(hist.at<float>(i))),
			Scalar(0, 0, 0), 1, 8, 0);
	}


	float factor = static_cast<float>(largest_hist_j * bin_w) / static_cast<float>(histImage.cols - 1);
	cout << "max value:  " << max_length << endl;
	cout << "peak value: " << max_length*factor << endl;


	circle(histImage, Point(largest_hist_j * bin_w, 0), 2, Scalar(255, 127, 0), 2);

	
	imshow("calcHist Demo", histImage);
	waitKey();
}


// TODO: fix camera bug where portrait mode crashes.
void take_screenshot(size_t num_cams_wide, const char *filename, const bool reverse_rows = false)
{
	get_points(150);

	// Set up Targa TGA image data.
	unsigned char  idlength = 0;
	unsigned char  colourmaptype = 0;
	unsigned char  datatypecode = 2;
	unsigned short int colourmaporigin = 0;
	unsigned short int colourmaplength = 0;
	unsigned char  colourmapdepth = 0;
	unsigned short int x_origin = 0;
	unsigned short int y_origin = 0;
	unsigned short int px = win_x*num_cams_wide;
	unsigned short int py = win_y*num_cams_wide;
	unsigned char  bitsperpixel = 24;
	unsigned char  imagedescriptor = 0;
	vector<char> idstring;

	size_t num_bytes = 3*px*py;
	vector<unsigned char> pixel_data(num_bytes);

	// Adjust some parameters for large screen format.
	bool temp_draw_control_list = draw_control_list;
	draw_control_list = false;

	float temp_outline_width = outline_width;
	outline_width = 12;

	vector<unsigned char> fbpixels(3*win_x*win_y);

	// Loop through subcameras.
	for(size_t cam_num_x = 0; cam_num_x < num_cams_wide; cam_num_x++)
	{
		for(size_t cam_num_y = 0; cam_num_y < num_cams_wide; cam_num_y++)
		{
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
		}
	}

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

	get_points(10);
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

	glClearColor(background_colour.x, background_colour.y, background_colour.z, 1);
	glClearDepth(1.0f);

	main_camera.Set(0, 0, camera_w, camera_fov, win_x, win_y, camera_near, camera_far);

	get_points(10);

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
		glTranslatef(0, -win_y, 0); // Neat. :)
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

	glFlush();
	glutSwapBuffers();
}

void keyboard_func(unsigned char key, int x, int y)
{
	switch(tolower(key))
	{
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
	case 'l':
		{
			take_screenshot(12, "screenshot.tga");
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
	if(false == disable_colouring)
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

	static const float rad_to_deg = 180.0f/pi;

	glPushMatrix();

	glTranslatef(camera_x_transform, camera_y_transform, 0);

	if (false == disable_colouring)
		glMaterialfv(GL_FRONT, GL_DIFFUSE, sphere_colour);

	for (size_t i = 0; i < 1/*pos.size()*/; i++)
	{
		for (size_t j = 0; j < pos[i].size() - 1; j++)
		{
			double t = j / static_cast<double>(pos[i].size() - 1);

			RGB rgb = HSBtoRGB(300 * t, 75, 100);

			float colour[] = { rgb.r / 255.0, rgb.g / 255.0, rgb.b / 255.0, 1.0f};

			glMaterialfv(GL_FRONT, GL_DIFFUSE, colour);

			static const float rad_to_deg = 180.0f / pi;

			vector_4 line = pos[i][j + 1] - pos[i][j];
			
			glPushMatrix();
			glTranslatef(pos[i][j].x, pos[i][j].y, pos[i][j].z);

			float line_len = line.length();
			line.normalize();
			
			float yaw = 0.0f;

			if (fabsf(line.x) < 0.00001 && fabsf(line.z) < 0.00001)
				yaw = 0.0f;
			else
				yaw = atan2f(line.x, line.z);

			float pitch = -atan2f(line.y, sqrt(line.x*line.x + line.z*line.z));

			glRotatef(yaw*rad_to_deg, 0.0f, 1.0f, 0.0f);
			glRotatef(pitch*rad_to_deg, 1.0f, 0.0f, 0.0f);

			gluCylinder(glu_obj, 0.005, 0.005, line_len, 20, 2);

			glPopMatrix();
		}
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

		glColor4f(0, 0, 0, 0.5);

		//		glColor4f(1, 0, 0, 0.25);
		glVertex3f(0, 0, 0);
		glVertex3f(1, 0, 0);
		//		glColor4f(0, 1, 0, 0.25);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 1, 0);
		//		glColor4f(0, 0, 1, 0.25);
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