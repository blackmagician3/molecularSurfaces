#ifndef CAMERA_HPP
#define CAMERA_HPP

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <vector>

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

// Define possible options for camera movement
enum Camera_Movement
{
	FORWARD,
	BACKWARD,
	LEFT,
	RIGHT,
	UPWARD,
	DOWNWARD
};
// Default camera values
const float YAW = -90.0f;
const float PITCH = 0.0f;
const float SPEED = 10.0f;
const float SENSITIVITY = 0.1f;
const float ZOOM = 10.0f;

const glm::vec3 FOCUS = glm::vec3(0.0f, 0.0f, 0.0f);

// An abstract camera class that processes input and calculates the corresponding Euler Angles, Vectors and Matrices for use in OpenGL
class Camera
{
public:
	// camera Attributes
	glm::vec3 Position;
	glm::vec3 Front;
	glm::vec3 Up;
	glm::vec3 Right;
	glm::vec3 WorldUp;
	glm::vec3 Focus;

	glm::vec3 mouseAxisVertical;
	glm::vec3 mouseAxisHorizontal;

	// camera options
	float MovementSpeed;
	float MouseSensitivity;
	float Zoom;

	// constructor with vectors
	Camera(glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f)) : Front(glm::vec3(0.0f, 0.0f, -1.0f)), MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM), Focus(FOCUS)
	{
		Position = position;
		WorldUp = up;
		xOffset = 0.0f;
		yOffset = 0.0f;

		updateCameraVectors();
	}
	// constructor with scalar values
	Camera(float posX, float posY, float posZ, float upX = 0.0f, float upY = 1.0f, float upZ = 0.0f) : Front(glm::vec3(0.0f, 0.0f, -1.0f)), MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM), Focus(Focus)
	{
		Position = glm::vec3(posX, posY, posZ);
		WorldUp = glm::vec3(upX, upY, upZ);

		updateCameraVectors();
	}

	// processes input received from any keyboard-like input system.
	void ProcessCameraMovement(Camera_Movement direction, float deltaTime)
	{
		float velocity = MovementSpeed * deltaTime;
		if (direction == FORWARD)
			Position += Front * velocity;
		if (direction == BACKWARD)
			Position -= Front * velocity;
		if (direction == LEFT)
		{
			Position -= Right * velocity;
			Focus -= Right * velocity;
		}
		if (direction == RIGHT)
		{
			Position += Right * velocity;
			Focus += Right * velocity;
		}
		if (direction == UPWARD)
		{
			Position += Up * velocity;
			Focus += Up * velocity;
		}
		if (direction == DOWNWARD)
		{
			Position -= Up * velocity;
			Focus -= Up * velocity;
		}
	}

	// processes input received from a mouse input system. Expects the offset value in both the x and y direction.
	void ProcessCameraRotation(float xoffset, float yoffset)
	{
		// TODO: optimize: perform overflow every 1000 or so rotations
		xOffset = xOffset + xoffset * MouseSensitivity;
		xOffset = overflowInRange(xOffset, 0, 2 * M_PI);
		yOffset = yOffset + yoffset * MouseSensitivity;
		yOffset = overflowInRange(yOffset, 0, 2 * M_PI);

		// express rotation with quaternions
		glm::quat rotateVertical = glm::angleAxis(glm::radians(-xoffset), mouseAxisVertical);
		glm::quat rotateHorizontal = glm::angleAxis(glm::radians(-yoffset), mouseAxisHorizontal);

		// cast
		glm::mat4 rot = glm::mat4(1.0);
		rot *= glm::mat4_cast(rotateVertical);
		rot *= glm::mat4_cast(rotateHorizontal);

		// apply rotation
		glm::vec4 newPos = rot * glm::vec4(Position.x, Position.y, Position.z, 0.0f);
		Position = glm::vec3(newPos.x, newPos.y, newPos.z);

		// update Front, Right and Up Vectors using the updated Euler angles
		updateCameraVectors();
	}
	float overflowInRange(float value, float min, float max)
	{

		float whole, frac, result;
		float range = max - min;
		frac = modf(value / range, &whole);
		result = frac * range + min;

		return result;
	}

	void intializeCameraPosition(float4 cameraStart)
	{
		Focus = glm::vec3(cameraStart.x, cameraStart.y, cameraStart.z);
		Position = glm::vec3(cameraStart.x, cameraStart.y, cameraStart.z + cameraStart.w);
		// Position = glm::vec3(-181.93, 113.00, -114.96);
		updateCameraVectors();
		updateMouseAxis();
	}

	void updateMouseAxis()
	{
		mouseAxisVertical = glm::normalize(Up);
		mouseAxisHorizontal = glm::normalize(Right);
	}

	// processes input received from a mouse scroll-wheel event. Only requires input on the vertical wheel-axis
	void ProcessMouseScroll(double yoffset)
	{
		yoffset *= 2 * MouseSensitivity;
		Zoom += (float)yoffset;
		if (Zoom < 1.0f)
			Zoom = 1.0f;
		if (Zoom > 90.0f)
			Zoom = 90.0f;
	}

	//////////////////////////////////
	///// getter (float4)
	//////////////////////////////////
	float4 getFocusAsFloat4()
	{
		float4 vector = make_float4(Focus.x, Focus.y, Focus.z, 0.0f);
		return vector;
	}
	float4 getCenterAsFloat4()
	{
		glm::vec3 center = Position + Front * Zoom;
		float4 vector = make_float4(center.x, center.y, center.z, 0.0f);
		return vector;
	}
	float4 getRightAsFloat4()
	{
		float4 vector = make_float4(Right.x, Right.y, Right.z, 0.0f);
		return vector;
	}
	float4 getUpAsFloat4()
	{
		float4 vector = make_float4(Up.x, Up.y, Up.z, 0.0f);
		return vector;
	}
	float4 getPosAsFloat4()
	{
		float4 vector = make_float4(Position.x, Position.y, Position.z, 0.0f);
		return vector;
	}
	float4 getFrontAsFloat4()
	{
		float4 vector = make_float4(Front.x, Front.y, Front.z, 0.0f);
		return vector;
	}

	void resetCamera()
	{
		// calculate remaining vectors
		Front = glm::normalize(Focus - Position);
		Right = glm::normalize(glm::cross(Front, WorldUp));
		Up = glm::normalize(glm::cross(Right, Front));
		Zoom = ZOOM;
		xOffset = 0.0f;
		yOffset = 0.0f;
	}

private:
	float xOffset;
	float yOffset;
	// calculates the front vector from the Camera's (updated) Euler Angles
	void updateCameraVectors()
	{
		// calculate the new Front vector
		// glm::vec3 front;
		// front.x = cos(glm::radians(Yaw)) * cos(glm::radians(Pitch));
		// front.y = sin(glm::radians(Pitch));
		// front.z = sin(glm::radians(Yaw)) * cos(glm::radians(Pitch));

		// calculate the new Front vector
		Front = glm::normalize(Focus - Position);

		// also re-calculate the Right and Up vector
		Right = glm::normalize(glm::cross(Front, WorldUp)); // normalize the vectors, because their length gets closer to 0 the more you look up or down which results in slower movement.
		Up = glm::normalize(glm::cross(Right, Front));
	}
};
#endif