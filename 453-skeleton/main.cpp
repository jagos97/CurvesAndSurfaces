#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <limits>
#include <functional>

#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Texture.h"
#include "Window.h"

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"


const float PI = 3.14159265359;

glm::vec3 red(1.0f, 0.f, 0.f);
glm::vec3 green(.0f, 1.f, 0.f);
glm::vec3 blue(.0f, 0.f, 1.f);
glm::vec3 yellow(1.0f, 1.f, 0.f);
glm::vec3 magenta(1.0f, 0.f, 1.f);
glm::vec3 cyan(.0f, 1.f, 1.f);
glm::vec3 black(0.0f, 0.f, 0.f);

int selected = -1;
int scene = 0;
int scene2 = 0;
int surfaceToUse = 0;

bool noneSelected;
bool someSelected;
bool notErased;
bool is3D = false;
float x = 0;
float y = 0;
float z = 2;
float yaw = PI / 2;
float pitch = 0;

float currentTime;
float deltaTime;
unsigned int* p3D;

int size;


// We gave this code in one of the tutorials, so leaving it here too
void updateGPUGeometry(GPU_Geometry &gpuGeom, CPU_Geometry const &cpuGeom) {
	gpuGeom.bind();
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setCols(cpuGeom.cols);
}

glm::vec3 Multiply(float x, glm::vec3 vec) {
	return glm::vec3(vec.x * x, vec.y * x, vec.z * x);
}

bool Close(glm::vec3 pos1, glm::vec2 pos2) {
	float x = pos2.x - pos1.x;
	float y = pos2.y - pos1.y;
	return sqrt(x * x + y * y) < 0.05f;
}

void AddNewPoint(CPU_Geometry &points, glm::vec2 pos) {
	points.verts.push_back(glm::vec3(pos, 0.f));
	points.cols.push_back(blue);
}

void MakeSelectedPointBlue(CPU_Geometry& controlPoints) {
	controlPoints.cols.clear();
	controlPoints.cols.resize(controlPoints.verts.size(), red);
	if (selected >= 0 && selected < controlPoints.verts.size()) {
		controlPoints.cols.at(selected) = blue;
	}
}

void DeleteSelectedPoint(CPU_Geometry& controlPoints) {
	controlPoints.verts.erase(controlPoints.verts.begin() + selected);
	controlPoints.cols.erase(controlPoints.cols.begin() + selected);
}


glm::mat4 MakeTranslationMatrixXY(float x, float y) {
	glm::mat4 translation(
		1.f, 0.f, 0.f, 0.f,
		0.f, 1.f, 0.f, 0.f,
		0.f, 0.f, 1.f, 0.f,
		x, y, 0.f, 1.f
	);
	return translation;
}

glm::mat4 MakeTranslationMatrixXYZ(float x, float y, float z) {
	glm::mat4 translation(
		1.f, 0.f, 0.f, 0.f,
		0.f, 1.f, 0.f, 0.f,
		0.f, 0.f, 1.f, 0.f,
		x, y, z, 1.f
	);
	return translation;
}



// EXAMPLE CALLBACKS
class Assignment3 : public CallbackInterface {

public:
	Assignment3()
	{
	}
	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if (!is3D) {
			if (action == GLFW_PRESS) {
				if (key == GLFW_KEY_X) {
					deletePoint = true;
				}
				else if (key == GLFW_KEY_RIGHT) {
					scene = 1;
				}
				else if (key == GLFW_KEY_LEFT) {
					scene = 0;
				}
				else if (key == GLFW_KEY_R) {
					reset = true;
				}
				else if (key == GLFW_KEY_SPACE) {
					if (!is3D) {
						make3D = true;
					}
				}
			}
		}
		else{
			if (action == GLFW_PRESS || action == GLFW_REPEAT) {
				if (key == GLFW_KEY_W) {
					glm::vec3 forward(glm::normalize(glm::vec3(cosf(yaw) * cosf(pitch), sinf(pitch), sinf(yaw) * cosf(pitch))));
					x -= rate * deltaTime * forward.x;
					y -= rate * deltaTime * forward.y;
					z -= rate * deltaTime * forward.z;
				}
				else if (key == GLFW_KEY_S) {
					glm::vec3 forward(glm::normalize(glm::vec3(cosf(yaw) * cosf(pitch), sinf(pitch), sinf(yaw) * cosf(pitch))));
					x += rate * deltaTime * forward.x;
					y += rate * deltaTime * forward.y;
					z += rate * deltaTime * forward.z;
				}
				else if (key == GLFW_KEY_A) {
					x += rate * deltaTime * cosf(yaw + PI / 2);
					z += rate * deltaTime * sinf(yaw+ PI/2);
				}
				else if (key == GLFW_KEY_D) {
					x -= rate * deltaTime * cosf(yaw + PI / 2);
					z -= rate * deltaTime * sinf(yaw + PI / 2);
				}
				else if (key == GLFW_KEY_X) {
					glm::vec3 forward(glm::normalize(glm::vec3(cosf(yaw) * cosf(pitch + PI / 2), sinf(pitch + PI / 2), sinf(yaw) * cosf(pitch + PI / 2))));
					x += rate * deltaTime * forward.x;
					y += rate * deltaTime * forward.y;
					z += rate * deltaTime * forward.z;
				}
				else if (key == GLFW_KEY_C) {
					glm::vec3 forward(glm::normalize(glm::vec3(cosf(yaw) * cosf(pitch + PI / 2), sinf(pitch + PI / 2), sinf(yaw) * cosf(pitch + PI / 2))));
					x -= rate * deltaTime * forward.x;
					y -= rate * deltaTime * forward.y;
					z -= rate * deltaTime * forward.z;
				}
				else if (key == GLFW_KEY_SPACE) {
					is3D = false;
				}
				else if (key == GLFW_KEY_RIGHT) {
					scene2++;
					if (scene2 > 3) {
						scene2 = 3;
					}
				}
				else if (key == GLFW_KEY_LEFT) {
					scene2--;
					if (scene2 < 0) {
						scene2 = 0;
					}
				}
				else if (key == GLFW_KEY_UP) {
					surfaceToUse++;
					if (surfaceToUse > 2) {
						surfaceToUse = 2;
					}
				}
				else if (key == GLFW_KEY_DOWN) {
					surfaceToUse --;
					if (surfaceToUse < 0) {
						surfaceToUse = 0;
					}
				}
			}
		}
		if (action == GLFW_PRESS && key == GLFW_KEY_P) {
			if (isWire) {
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			}
			else glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			isWire = !isWire;
		}
	}
	virtual void mouseButtonCallback(int button, int action, int mods) {
		if (button == GLFW_MOUSE_BUTTON_LEFT) {
			if (action == GLFW_PRESS) {
				GLmouse();
				leftPressed = true;
				firstClick = true;
				needToSelect = true;
				notErased = true;
			}
			else if (action == GLFW_RELEASE) {
				leftPressed = false;
			}
		}	
	}
	virtual void cursorPosCallback(double xpos, double ypos) {
		if (leftPressed && is3D) {
			float deltax, deltay;
			if (firstClick) {
				firstClick = false;
			}
			else {
				deltax = xpos - mousePos.x;
				deltay = mousePos.y - ypos;
				yaw += deltax * scale;
				pitch -= deltay * scale;
				if (pitch >= PI / 2) {
						pitch = PI / 2 - 0.02;
				}
				else if (pitch < -PI / 2) pitch = -PI / 2 + 0.02;
			}
		}
		mousePos.x = xpos;
		mousePos.y = ypos;
	}
	virtual void scrollCallback(double xoffset, double yoffset) {
	}
	virtual void windowSizeCallback(int width, int height) {
		sWidth = width / 2;
		sHeight = height / 2;
		
		// The CallbackInterface::windowSizeCallback will call glViewport for us
		CallbackInterface::windowSizeCallback(width,  height);
	}
	void GLmouse() {
		clickPos = glm::vec2(mousePos.x / sWidth - 1, mousePos.y / -sHeight + 1);
	}
	bool GetLeftPressed() {
		return leftPressed;
	}
	bool ShouldDeletePoint() {
		return deletePoint;
	}
	void PointDeleted() {
		deletePoint = false;
	}

	bool NeedToSelect() {
		return needToSelect;
	}
	void SetSelected() {
		needToSelect = false;
	}
	bool ShouldReset() {
		return reset;
	}
	void Resetted() {
		reset = false;
	}
	bool Make3D() {
		return make3D;
	}
	void Made3D() {
		make3D = false;
	}
	glm::vec2 GetClickPos() {
		return clickPos;
	}
private:
	float scale = 0.01;
	float const rate = 80;
	bool firstClick;
	glm::vec2 mousePos = glm::vec2(1.0f);
	glm::vec2 clickPos;
	bool leftPressed = false;
	bool needToSelect = false;
	bool deletePoint = false;
	bool reset = false;
	bool make3D = false;
	int sWidth = 400;
	int sHeight = 400;
	bool isWire = false;
};

void AddToVector(std::vector<glm::vec3> &toadd, const std::vector<glm::vec3> adding) {
	for (int i = 0; i < adding.size(); i++) {
		toadd.push_back(adding.at(i));
	}
}

void MakeSame(CPU_Geometry &a, std::vector<glm::vec3> b) {
	a.verts.clear();
	a.verts.resize(b.size());
	for (int i = 0; i < b.size(); i++) {
		a.verts.at(i) = glm::vec3(b.at(i).x, b.at(i).y, b.at(i).z);
	}


}



void MakeBezierCurve(CPU_Geometry cp, CPU_Geometry &line) {
	int d = cp.verts.size() - 1;
	line.verts.clear();
	CPU_Geometry temp;
	MakeSame(temp, cp.verts);
	for (float u = 0; u <= 1; u = u + 0.02) {
		for (int i = 1; i <= d; i++) {
			for (int j = 0; j <= d - i; j++) {
				temp.verts.at(j) = Multiply(1-u,temp.verts.at(j)) + Multiply(u, temp.verts.at(j + 1));
			}
		}
		if(!cp.verts.empty())line.verts.push_back(temp.verts.at(0));
		MakeSame(temp, cp.verts);
	}

}

void MakeBezierCurveWithVectors( std::vector<glm::vec3> cp, std::vector<glm::vec3> &line) {
	int d = cp.size() - 1;
	line.clear();
	CPU_Geometry temp;
	MakeSame(temp, cp);
	for (float u = 0; u <= 1; u = u + 0.02) {
		for (int i = 1; i <= d; i++) {
			for (int j = 0; j <= d - i; j++) {
				temp.verts.at(j) = Multiply(1 - u, temp.verts.at(j)) + Multiply(u, temp.verts.at(j + 1));
			}
		}
		if (!cp.empty())line.push_back(temp.verts.at(0));
		MakeSame(temp, cp);
	}
}


void MakeBSpline(CPU_Geometry const cp, CPU_Geometry& line) {
	int j = 0;
	int last = cp.verts.size() - 1;
	line.verts.clear();
	CPU_Geometry temp;
	MakeSame(temp, cp.verts);
	for (int d = 0; d < 5; d++) {
		if (last >= 1) {
			line.verts.push_back(temp.verts.at(0));
			line.verts.push_back(Multiply(0.5f, temp.verts.at(0)) + Multiply(0.5f, temp.verts.at(1)));
			for (int i = 1; i < temp.verts.size() - 2; i++) {
				line.verts.push_back(Multiply(0.75f, temp.verts.at(i)) + Multiply(0.25f, temp.verts.at(i + 1)));
				line.verts.push_back(Multiply(0.25f, temp.verts.at(i)) + Multiply(0.75f, temp.verts.at(i + 1)));
				j += 2;
			}
			//line.verts.push_back(Multiply(0.75f, cp.verts.at(last)) + Multiply(0.25f, cp.verts.at(0)));
			line.verts.push_back(Multiply(0.5f, temp.verts.at(last)) + Multiply(0.5f, temp.verts.at(last - 1)));
			line.verts.push_back(temp.verts.at(last));
			if (d < 4) {
				MakeSame(temp, line.verts);
				line.verts.clear();
			}
		}
		last = temp.verts.size()-1;
	}

	

}


void MakeBSplineWithVectors(std::vector<glm::vec3> cp, std::vector<glm::vec3>& line) {
	int j = 0;
	int last = cp.size() - 1;
	line.clear();
	CPU_Geometry temp;
	MakeSame(temp, cp);
	for (int d = 0; d < 5; d++) {
		if (last >= 1) {
			line.push_back(temp.verts.at(0));
			line.push_back(Multiply(0.5f, temp.verts.at(0)) + Multiply(0.5f, temp.verts.at(1)));
			for (int i = 1; i < temp.verts.size() - 2; i++) {
				line.push_back(Multiply(0.75f, temp.verts.at(i)) + Multiply(0.25f, temp.verts.at(i + 1)));
				line.push_back(Multiply(0.25f, temp.verts.at(i)) + Multiply(0.75f, temp.verts.at(i + 1)));
				j += 2;
			}
			line.push_back(Multiply(0.5f, temp.verts.at(last)) + Multiply(0.5f, temp.verts.at(last - 1)));
			line.push_back(temp.verts.at(last));
			if (d < 4) {
				MakeSame(temp, line);
				line.clear();
			}
		}
		last = temp.verts.size() - 1;
	}
	
}



void MakeCurve(CPU_Geometry cp, CPU_Geometry& line) {
	if (scene == 0)  MakeBezierCurve(cp, line);
	else MakeBSpline(cp, line);
}


void PutPointsIntoSurface(CPU_Geometry &surface, std::vector<std::vector<glm::vec3>>& const points) {
	static unsigned int temp[500000];
	int index = 0;
	int last = points.size() - 1;
	int last2 = points.at(0).size() - 1;
	int sizeOfx = points.size();
	int sizeOfy = points.at(0).size();

	for (auto i = points.begin(); i < points.end() -1 ; i++) {
		for (int j = 0; j < (*i).size() - 1; j++) {
			temp[index] = (i - points.begin()) * sizeOfy + j +1;
			temp[index + 1] = (i - points.begin()) * sizeOfy + j;
			temp[index + 2] = (i - points.begin() + 1) * sizeOfy + j + 1;
			temp[index + 3] = (i - points.begin() + 1) * sizeOfy + j + 1;
			temp[index + 4] = (i - points.begin()) * sizeOfy + j;
			temp[index + 5] = (i - points.begin() + 1) * sizeOfy + j;
			index += 6;



		}
	}

	if (scene2 == 1) {
		for (int j = 0; j < last2; j++) {
			temp[index] = last * sizeOfy + j + 1;
			temp[index + 1] = last * sizeOfy + j;
			temp[index + 2] = j + 1;
			temp[index + 3] = j + 1;
			temp[index + 4] = last * sizeOfy + j;
			temp[index + 5] = j;
			index += 6;
			
		}
	}

	//std::cout << index << std::endl;			//comment me out
	size = index;
	p3D = temp;


}

void MakeSurface(CPU_Geometry lines, CPU_Geometry& surface) {
	surface.verts.clear();
	if (lines.verts.size() > 0) {
		int m = 0;
		std::vector<std::vector<glm::vec3>> points;
		for (float u = 0; u <= 2 * PI; u += 0.1) {
			points.push_back(std::vector<glm::vec3>());
			for (auto i = lines.verts.begin(); i < lines.verts.end(); i++) {
				points.at(m).push_back(glm::vec3((*i).x * cosf(u), (*i).y, (*i).x * sinf(u)));
				surface.verts.push_back(glm::vec3((*i).x * cosf(u), (*i).y, (*i).x * sinf(u))   );
			}
			m++;
		}
		PutPointsIntoSurface(surface, points);
		surface.cols.clear();
		surface.cols.resize(surface.verts.size(), black);
	}
}

void MakeBezierSurface(std::vector<std::vector<glm::vec3>> & cp, CPU_Geometry &surface) {
	std::vector<std::vector<glm::vec3>> R;
	std::vector<std::vector<glm::vec3>> Q;
	for (int i = 0; i < cp.size(); i++) {
		R.push_back(std::vector<glm::vec3>());
		MakeBezierCurveWithVectors(cp.at(i), R.at(i));
	}

	std::vector<std::vector<glm::vec3>> temp;
	for (int j = 0; j < R.at(0).size(); j++ ) {
		temp.push_back(std::vector<glm::vec3>());
		for (int i = 0; i < R.size(); i++) {
			temp.at(j).push_back(R.at(i).at(j));
		}
	}
	for (int i = 0; i < temp.size(); i++) {
		Q.push_back(std::vector<glm::vec3>());
		MakeBezierCurveWithVectors(temp.at(i), Q.at(i));
		AddToVector(surface.verts, Q.at(i));

	}
	PutPointsIntoSurface(surface, Q);
	surface.cols.clear();
	surface.cols.resize(surface.verts.size(), black);
}

void MakeChaikinSurface(std::vector<std::vector<glm::vec3>> &cp, CPU_Geometry &surface) {
	std::vector<std::vector<glm::vec3>> R;
	std::vector<std::vector<glm::vec3>> Q;
	for (int i = 0; i < cp.size(); i++) {
		R.push_back(std::vector<glm::vec3>());
		MakeBSplineWithVectors(cp.at(i), R.at(i));
	}
	std::vector<std::vector<glm::vec3>> temp;
	for (int j = 0; j < R.at(0).size(); j++) {
		temp.push_back(std::vector<glm::vec3>());
		for (int i = 0; i < R.size(); i++) {
			temp.at(j).push_back(R.at(i).at(j));
		}
	}
	for (int i = 0; i < temp.size(); i++) {
		Q.push_back(std::vector<glm::vec3>());
		MakeBSplineWithVectors(temp.at(i), Q.at(i));
		AddToVector(surface.verts, Q.at(i) );
	}
	PutPointsIntoSurface(surface, Q);
	surface.cols.clear();
	surface.cols.resize(surface.verts.size(), black);

}

void MakeTensor(std::vector<std::vector<glm::vec3>> &cp, CPU_Geometry &surface) {
	surface.verts.clear();
	if (scene2 == 2) MakeBezierSurface(cp, surface);
	else MakeChaikinSurface(cp, surface);

	

}


void UpdateGPUS(GPU_Geometry &pointsGPUGeom, CPU_Geometry &controlPoints, GPU_Geometry &cpLines, GPU_Geometry &linesGPUGeom, CPU_Geometry &lines) {
	updateGPUGeometry(pointsGPUGeom, controlPoints);
	controlPoints.cols.clear();
	controlPoints.cols.resize(controlPoints.verts.size(), green);
	updateGPUGeometry(cpLines, controlPoints);
	updateGPUGeometry(linesGPUGeom, lines);
}



/*Can handle any bezier curve but only bpline with up to 43 control points
*/
unsigned int* MakeIndicesForLines() {
	static unsigned int r[1316];


	for (int i = 0; i < 1316; i++) {
		r[i] = i;
	}
	return r;
}

unsigned int* MakeIndicesForSurfaceLines() {
	static unsigned int r[24] = { 0,1, 0,3,  1,2, 1,4, 2,5, 5,4  ,4,3  ,3,6, 6,7,  7,4, 7,8, 8,5};

	return r;

}


int main() {
	Log::debug("Starting main");

	// WINDOW
	glfwInit();
	Window window(800, 800, "CPSC 453"); // can set callbacks at construction if desired


	GLDebug::enable();

	auto a3 = std::make_shared<Assignment3>();
	window.setCallbacks(a3);
	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");

	// The current CPU_Geometry and GPU_Geometry classes are defined in
	// Geometry.h/Geometry.cpp They will work for this assignment, but for some of
	// the bonuses you may have to modify them.
	CPU_Geometry controlPoints;
	CPU_Geometry quads;
	controlPoints.verts.push_back(glm::vec4{ 0, 0.5, 0,1 });
	controlPoints.verts.push_back(glm::vec4{ 0.5, -0.5, 0,1 });

	glm::mat4 identity = glm::mat4(1.f);

	glm::vec3 camPos(x, y, z);
	glm::vec3 center(0, 0, z - 1);
	glm::vec3 up(0, 1, 0);

	glm::mat4 view(1.f);

	a3->GLmouse();
	view = glm::lookAt(camPos, center, up);


	float fov = PI / 2;
	float nearPane = 0.01f;
	float farPane = 10.f;

	glm::mat4 per = glm::perspective(fov,
		1.0f,
		nearPane,
		farPane
	);


	GPU_Geometry pointsGPUGeom;
	controlPoints.cols.clear();
	controlPoints.cols.resize(controlPoints.verts.size(), red);
	updateGPUGeometry(pointsGPUGeom, controlPoints);


	GPU_Geometry cpLines;
	controlPoints.cols.clear();
	controlPoints.cols.resize(controlPoints.verts.size(), green);
	updateGPUGeometry(cpLines, controlPoints);

	CPU_Geometry surface;
	GPU_Geometry surfaceGeom;






	CPU_Geometry lines;
	MakeCurve(controlPoints, lines);
	lines.cols.resize(lines.verts.size(), glm::vec3{ 0.0, 0.0, 0.0 });
	GPU_Geometry linesGPUGeom;
	updateGPUGeometry(linesGPUGeom, lines);
	glPointSize(10.0f);

	// Note this call only work on some systems, unfortunately.
	// In order for them to work, you have to comment out line 60
	// If you're on a mac, you can't comment out line 60, so you
	// these will have no effect. :(
	// glLineWidth(5.0f);


	int currentScene = scene;
	int currentScene2 = scene2;
	int currentSurface = surfaceToUse;
	


	glm::vec3 s11(-1.f, 0.f, -2.f);
	glm::vec3 s12(-1.f, 2.f, -1.f);
	glm::vec3 s13(-1.f, 0.f, 0.f);
	glm::vec3 s21(0.f, 0.f, -2.f);
	glm::vec3 s22(0.f, -2.f, -1.f);
	glm::vec3 s23(1.f, 0.f, 0.f);
	glm::vec3 s31(1.f, 0.f, -2.f);
	glm::vec3 s32(1.f, 1.f, -0.5f);
	glm::vec3 s33(2.f, 2.f, 0.5f);
	std::vector<glm::vec3> s1 = { s11,s12,s13 };
	std::vector<glm::vec3> s2 = { s21,s22,s23 };
	std::vector<glm::vec3> s3 = { s31,s32,s33 };
	std::vector<glm::vec3> s4 = { s11, s12, s13, s21,s22,s23,s31,s32,s33 };

	std::vector<std::vector<glm::vec3>> surfacePoints2 = { s1,s2,s3 };

	glm::vec3 a11(-2.f, 2.f, -2.f);
	glm::vec3 a12(-1.f, 0.f, -1.f);
	glm::vec3 a13(-2.f, 2.f, 0.f);
	glm::vec3 a21(0.f, 0.f, -1.f);
	glm::vec3 a22(0.f, 5.f, -.5f);
	glm::vec3 a23(2.f, 0.f, 0.f);
	glm::vec3 a31(1.f, 1.f, -2.f);
	glm::vec3 a32(2.f, 2.f, -0.5f);
	glm::vec3 a33(3.f, 0.f, 0.5f);
	std::vector<glm::vec3> a1 = { a11,a12,a13 };
	std::vector<glm::vec3> a2 = { a21,a22,a23 };
	std::vector<glm::vec3> a4 = { a31,a32,a33 };
	std::vector<glm::vec3> a5 = { a11, a12, a13, a21,a22,a23,a31,a32,a33};

	std::vector<std::vector<glm::vec3>> surfacePoints = { a1,a2,a4 };


	glm::vec3 b11(-1.f, 2.f, -2.f);
	glm::vec3 b12(-1.f, 0.f, 0.f);
	glm::vec3 b13(2.f, 2.f, 2.f);
	glm::vec3 b21(0.f, 0.f, -2.f);
	glm::vec3 b22(0.f, 1.f, -1.f);
	glm::vec3 b23(2.f, 0.f, 0.f);
	glm::vec3 b31(1.f, 1.f, -2.f);
	glm::vec3 b32(2.f, 2.f, -0.5f);
	glm::vec3 b33(3.f, 0.f, 0.5f);
	std::vector<glm::vec3> b1 = { b11,b12,b13 };
	std::vector<glm::vec3> b2 = { b21,b22,b23 };
	std::vector<glm::vec3> b3 = { b31,b32,b33 };
	std::vector<glm::vec3> b4 = { b11, b12, b13, b21,b22,b23,b31,b32,b33 };

	std::vector<std::vector<glm::vec3>> surfacePoints3 = { b1,b2,b3 };

	std::vector<std::vector<glm::vec3>> surfacesCP[3] = { surfacePoints, surfacePoints2, surfacePoints3 };
	std::vector<glm::vec3> surfaceLines[3] = { a5,s4, b4 };

	CPU_Geometry linesToDrawForSurface;
	GPU_Geometry GPULinesToDraw;

	unsigned int* pLines = MakeIndicesForLines();
	unsigned int* pLines2 = MakeIndicesForSurfaceLines();

	glEnable(GL_DEPTH_TEST);
	currentTime = glfwGetTime();
	// RENDER LOOP
	while (!window.shouldClose()) {
		deltaTime = glfwGetTime() - currentTime;
		currentTime = glfwGetTime();

		glfwPollEvents();




		

		GLint viewLoc = glGetUniformLocation(shader, "view");
		GLint perLoc = glGetUniformLocation(shader, "per");

		
		//If the user is clicking
		if (a3->GetLeftPressed() && notErased && !is3D) {
			//If no point has been selected
			if (a3->NeedToSelect()) {
				a3->GLmouse();
				noneSelected = true;
				//Loop to check if any of the control points have been selected
				for (int i = 0; i < controlPoints.verts.size(); i++) {
					//If clicked on a control point then select that one (if 2 or more are considered "Close", it will pick the earlier one)
					if (Close(controlPoints.verts.at(i), a3->GetClickPos())) {
						selected = i;
						noneSelected = false;
						break;
					}
				}
				//If user didn't click on an existing point then create a new one
				if (noneSelected) {
					AddNewPoint(controlPoints, a3->GetClickPos());
					MakeCurve(controlPoints, lines);
					selected = controlPoints.verts.size() - 1;
				}
				a3->SetSelected();
			}
			//If point has been selected
			else {
				a3->GLmouse();
				//Move selected point if it exists (Hasn't been deleted or reset)
				if (!controlPoints.verts.empty() && selected != -1) controlPoints.verts.at(selected) = glm::vec3(a3->GetClickPos().x, a3->GetClickPos().y,0.f);
				MakeCurve(controlPoints, lines);
			}
			MakeSelectedPointBlue(controlPoints);
			UpdateGPUS(pointsGPUGeom, controlPoints, cpLines, linesGPUGeom, lines);
			//std::cout << lines.verts.size() << " with cp size " << controlPoints.verts.size()<< std::endl;							///comment me out!!!!!!
		}
		//Delete a point
		if (a3->ShouldDeletePoint() && !is3D) {
			//Delete the point if point exists
			if (selected != -1) {
				DeleteSelectedPoint(controlPoints);
				controlPoints.cols.clear();
				controlPoints.cols.resize(controlPoints.verts.size(), red);
				selected = -1;
			}
			//make new curve
			if (controlPoints.verts.size() != 0) {
				MakeCurve(controlPoints, lines);
			}
			a3->PointDeleted();
			notErased = false;
			MakeSelectedPointBlue(controlPoints);
			UpdateGPUS(pointsGPUGeom, controlPoints, cpLines, linesGPUGeom, lines);
		}

		//change scene
		if ((currentScene != scene || currentScene2 != scene2) ) {
			currentScene = scene;
			if (scene == 0) MakeBezierCurve(controlPoints, lines);
			else if (scene == 1) MakeBSpline(controlPoints, lines);
			MakeSelectedPointBlue(controlPoints);
			UpdateGPUS(pointsGPUGeom, controlPoints, cpLines, linesGPUGeom, lines);
			if (is3D) {
				currentScene2 = scene2;
				if (scene2 == 1) {
					MakeSurface(lines, surface);
					updateGPUGeometry(surfaceGeom, surface);

				}
				else if (scene2 > 1){
					linesToDrawForSurface.verts.clear();
					linesToDrawForSurface.verts = surfaceLines[surfaceToUse];
					linesToDrawForSurface.cols.clear();
					linesToDrawForSurface.cols.resize(linesToDrawForSurface.verts.size(), blue);
					MakeTensor(surfacesCP[surfaceToUse], surface);
					updateGPUGeometry(surfaceGeom, surface);
				}
			}
		}
		//reset scene
		if (a3->ShouldReset() && !is3D) {
			controlPoints.verts.clear();
			lines.verts.clear();
			UpdateGPUS(pointsGPUGeom, controlPoints, cpLines, linesGPUGeom, lines);
			a3->Resetted();
			selected = -1;
		}

		if (a3->Make3D()) {

			if (currentScene2 == 1) {
				MakeSurface(lines, surface);
			}
			else if(scene2 >1) {
				linesToDrawForSurface.verts.clear();
				linesToDrawForSurface.verts = surfaceLines[surfaceToUse];
				linesToDrawForSurface.cols.clear();
				linesToDrawForSurface.cols.resize(linesToDrawForSurface.verts.size(), blue);
				updateGPUGeometry(GPULinesToDraw, linesToDrawForSurface);
				MakeTensor(surfacesCP[surfaceToUse], surface);
			}
			updateGPUGeometry(surfaceGeom, surface);
			is3D = true;
			a3->Made3D();
		}

		
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_FRAMEBUFFER_SRGB);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		shader.use();


		if (is3D) {
			
			
			camPos = glm::vec3(x, y, z);
			glm::vec3 front = glm::normalize(glm::vec3(cosf(yaw) * cosf(pitch), sinf(pitch), sinf(yaw) * cosf(pitch)));
			view = glm::lookAt(camPos, camPos - front, up);

			glUniformMatrix4fv(viewLoc,
				1,
				false,
				glm::value_ptr(view)
			);

			glUniformMatrix4fv(perLoc,
				1,
				false,
				glm::value_ptr(per)
			);
		}
		else {
			glUniformMatrix4fv(viewLoc,
				1,
				false,
				glm::value_ptr(identity)
			);

			glUniformMatrix4fv(perLoc,
				1,
				false,
				glm::value_ptr(identity)
			);

		}

		if (currentSurface != surfaceToUse) {
			if (scene2 >= 2) {
				linesToDrawForSurface.verts.clear();
				linesToDrawForSurface.verts = surfaceLines[surfaceToUse];
				linesToDrawForSurface.cols.clear();
				linesToDrawForSurface.cols.resize(linesToDrawForSurface.verts.size(), blue);
				MakeTensor(surfacesCP[surfaceToUse], surface);
			}
			currentSurface = surfaceToUse;
		}
		
		if (!is3D || scene2 == 0) {
			pointsGPUGeom.bind();
			glDrawArrays(GL_POINTS, 0, GLsizei(controlPoints.verts.size()));


			linesGPUGeom.bind();
			glDrawElements(GL_LINE_STRIP, lines.verts.size(), GL_UNSIGNED_INT, pLines);
		}
		else {
			if (! (scene2 == 1 && controlPoints.verts.size() == 0)) {
				updateGPUGeometry(surfaceGeom, surface);
				glDrawElements(GL_TRIANGLES, size, GL_UNSIGNED_INT, p3D);
			}
			if (scene2 >= 2) {
				updateGPUGeometry(GPULinesToDraw, linesToDrawForSurface);
				glDrawElements(GL_LINES, 24,  GL_UNSIGNED_INT, pLines2 );

				updateGPUGeometry(GPULinesToDraw, linesToDrawForSurface);
				glDrawArrays(GL_POINTS, 0, GLsizei(linesToDrawForSurface.verts.size()));
			}

		}

		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		window.swapBuffers();
	}

	glfwTerminate();
	return 0;
}
