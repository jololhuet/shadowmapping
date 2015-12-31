#version 330 core
uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

uniform bool useCsm;

uniform mat4 depth_vp_offset;
uniform mat4 depth_vp_offset0;
uniform mat4 depth_vp_offset1;
uniform mat4 depth_vp_offset2;
uniform mat4 depth_vp_offset3;
uniform mat4 depth_vp_offset4;
uniform mat4 depth_vp_offset5;
uniform mat4 depth_vp_offset6;
uniform mat4 depth_vp_offset7;

uniform float farPlane[8];

uniform float texRatio;
uniform bool has_normal;

uniform bool useNormalOffset;

in vec3 vpoint;
in vec3 vnormal;
in vec2 vtexcoord;
out vec3 n;
out vec2 uv;
out vec4 shadow_coord;
out vec4 vpoint_MV;

out float distToCamera;

vec4 computeShadowCoord (vec4 shadow_pos) {
    vec4 sc;
    if (useCsm) {
        if (distToCamera <= farPlane[0]) {
            sc = depth_vp_offset0 * model * shadow_pos;
        } else if (distToCamera <= farPlane[1]) {
            sc = depth_vp_offset1 * model * shadow_pos;
        } else if (distToCamera <= farPlane[2]) {
            sc = depth_vp_offset2 * model * shadow_pos;
        } else if (distToCamera <= farPlane[3]) {
            sc = depth_vp_offset3 * model * shadow_pos;
        } else if (distToCamera <= farPlane[4]) {
            sc = depth_vp_offset4 * model * shadow_pos;
        } else if (distToCamera <= farPlane[5]) {
            sc = depth_vp_offset5 * model * shadow_pos;
        } else if (distToCamera <= farPlane[6]) {
            sc = depth_vp_offset6 * model * shadow_pos;
        } else if (distToCamera <= farPlane[7]) {
            sc = depth_vp_offset7 * model * shadow_pos;
        } else {
            vec3(0);
        }
    } else {
        sc = depth_vp_offset * model * shadow_pos;
    }

    return sc;
}

void main() {
    // Inverse transpose
    n = (inverse(transpose(model)) * vec4(vnormal,0.0)).xyz;

    vec4 shadow_pos = vec4(vpoint,1.0);

    vec4 vpoint_mv = view * model * vec4(vpoint, 1.0);
    vpoint_MV = vpoint_mv;

    gl_Position = projection * vpoint_mv;
    distToCamera = gl_Position.z;

    if (useNormalOffset) {
        float normalOffsetScale = 0.0015;
        vec3 normalOffset = normalize(n) * normalOffsetScale;

        shadow_pos += vec4(normalOffset, 0.0);
    }

    shadow_coord = computeShadowCoord(shadow_pos);

    uv = vtexcoord;
}
