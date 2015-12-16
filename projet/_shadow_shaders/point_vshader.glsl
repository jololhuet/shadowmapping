#version 330 core

uniform mat4 view;
uniform mat4 projection;
uniform mat4 vp;

uniform vec3 position;

void main(void)
{
    gl_PointSize = 6.0;
    gl_Position = vp * vec4(position, 1.0);
}

