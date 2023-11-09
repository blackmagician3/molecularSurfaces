#version 330 core

out vec4 FragColor;
in vec2 texCoord;

uniform sampler2D texture_1;
uniform sampler2D texture_2;

void main()
{
    // FragColor = mix(texture(texture_1, texCoord), texture(texture_2, texCoord), 0.0);
    FragColor = texture(texture_1, texCoord);
}