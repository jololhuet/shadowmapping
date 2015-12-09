#version 330 core
out vec3 color;

uniform int light_type;
uniform vec3 light_pos = vec3(0);

in vec4 shadow_coord;

void main() {
    color = vec3(shadow_coord.z / shadow_coord.w);
}
