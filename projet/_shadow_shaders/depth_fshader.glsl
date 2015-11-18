#version 330 core
out vec3 color;

uniform int light_type;
uniform vec3 light_pos = vec3(0);

in vec4 shadow_coord;

void main() {
    float distanceToLight = 1.0;

    color = vec3(shadow_coord.z / shadow_coord.w);
//    if (distanceToLight > 1) {
//        color = vec3(1.0, 0.0, 0.0);
//    }
//        color = vec3(distanceToLight);
}
