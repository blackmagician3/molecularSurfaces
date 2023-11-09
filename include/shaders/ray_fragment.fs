#version 330 core
out vec4 FragColor;


in vec3 ourColor;
in vec2 TexCoord;

//texture sampler
uniform sampler2D texture1;
uniform sampler2D texture2;

// color for opacity
vec4 aColor = vec4(1.0, 1.0, 1.0, 1.0);

void main()
{
	// insert code here
	vec2 Coord = vec2(TexCoord.xy);

	aColor = vec4(1.0, 0.0, 0.0, 1.0);
	FragColor = mix(texture(texture1, Coord), texture(texture2, Coord), 0.0)*aColor;
}