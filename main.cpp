// Shawn Halayka -- shalayka@gmail.com
// June 26, 2010
//
// This code and data is in the public domain.

// This code uses the OpenGL 1.1 and GLUT APIs.
// - For *nix: use the system package manager to install them, if need be.
// - For Windows: use Nate Robins' GLUT implementation, or freeglut, or OpenGLUT, etc.


#include "main.h"

int main(int argc, char **argv)
{	
	glutInit(&argc, argv);
	init_opengl(win_x, win_y);
	glutReshapeFunc(reshape_func);
	glutIdleFunc(idle_func);
	glutDisplayFunc(display_func);
	glutKeyboardFunc(keyboard_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutPassiveMotionFunc(passive_motion_func);
	//glutIgnoreKeyRepeat(1);

	glutMainLoop();

	glutDestroyWindow(win_id);

	return 0;
}
