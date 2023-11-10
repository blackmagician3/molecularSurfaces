#ifndef CAMERA_H
#define CAMERA_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <vector>

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
const float SPEED = 5.0f;
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

	// camera options
	float MovementSpeed;
	float MouseSensitivity;
	float Zoom;

	// constructor with vectors
	Camera(glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f)) : Front(glm::vec3(0.0f, 0.0f, -1.0f)), MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM), Focus(FOCUS)
	{
		Position = position;
		WorldUp = up;

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
	void ProcessCameraRotation(float xoffset, float yoffset, GLboolean constrainPitch = true)
	{
		xoffset *= MouseSensitivity;
		yoffset *= MouseSensitivity;

		// step1: find rotation axis

		glm::vec3 axisVerticalMouse = glm::normalize(Focus + Up);
		glm::vec3 axisHorizontalMouse = glm::normalize(Focus + Right);

		// express rotation with quaternions
		glm::quat rotateVertical = glm::angleAxis(glm::radians(-xoffset), axisVerticalMouse);
		glm::quat rotateHorizontal = glm::angleAxis(glm::radians(-yoffset), axisHorizontalMouse);

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

	// processes input received from a mouse scroll-wheel event. Only requires input on the vertical wheel-axis
	void ProcessMouseScroll(double yoffset)
	{
		yoffset *= MouseSensitivity;
		Zoom += (float)yoffset;
		if (Zoom < 1.0f)
			Zoom = 1.0f;
		if (Zoom > 45.0f)
			Zoom = 45.0f;
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

	void resetCamera(glm::vec3 cameralocation)
	{
		// set position to starting position
		// TODO: set starting position (z-component based on molecule dimensions)
		Position = cameralocation;

		// reset focus to the center of the world view
		Focus = FOCUS;

		// calculate remaining vectors
		Front = glm::normalize(Focus - Position);
		Right = glm::normalize(glm::cross(Front, WorldUp));
		Up = glm::normalize(glm::cross(Right, Front));
		Zoom = ZOOM;
	}

private:
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