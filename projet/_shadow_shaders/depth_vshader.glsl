#version 330 core
uniform mat4 view;
uniform mat4 model;

uniform int light_type;

uniform mat4 depth_vp_offset;

in vec3 vpoint;

out vec4 shadow_coord;

void main() {
    shadow_coord = depth_vp_offset * model * vec4(vpoint,1.0);

    gl_Position = depth_vp_offset * model * vec4(vpoint,1.0) -0.5;
}
