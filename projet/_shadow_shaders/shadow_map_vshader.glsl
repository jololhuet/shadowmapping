#version 330 core
uniform bool isVsm;
uniform mat4 depth_vp;
uniform mat4 depth_vp_offset;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

in vec3 vnormal;
in vec3 vpoint;

out vec4 vpoint_mv;

void main() {
    gl_Position = depth_vp * model * vec4(vpoint, 1.0);
    vpoint_mv = gl_Position;
}
