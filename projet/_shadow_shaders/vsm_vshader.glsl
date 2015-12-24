#version 330 core
uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;
uniform mat4 depth_vp_offset;

uniform sampler2D normalTex;
uniform float texRatio;
uniform bool has_normal;

in vec3 vpoint;
in vec3 vnormal;
in vec2 vtexcoord;
out vec3 n;
out vec2 uv;
out vec4 shadow_coord;
out vec4 vpoint_MV;
out float distToCamera;

void main() {

    n = (inverse(transpose(model)) * vec4(vnormal,0.0)).xyz;

    if (has_normal) {
        n = mat3(model) * texture2D(normalTex, (uv*texRatio)).rgb;
        n = vec3(n.g, n.r, n.b);
    }

    vec4 shadow_pos = vec4(vpoint,1.0);

    vec4 vpoint_mv = view * model * vec4(vpoint, 1.0);
    vpoint_MV = vpoint_mv;

    shadow_coord = depth_vp_offset * model * shadow_pos;

    gl_Position = projection * vpoint_mv;
    distToCamera = gl_Position.z;
    uv = vtexcoord;
}
