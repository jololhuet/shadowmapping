#version 330 core
uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;
uniform mat4 depth_vp_offset;

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

void main() {
    // Inverse transpose
    n = (inverse(transpose(model)) * vec4(vnormal,0.0)).xyz;

    vec4 shadow_pos = vec4(vpoint,1.0);

    if (useNormalOffset) {
        float normalOffsetScale = 0.0015;
        vec3 normalOffset = normalize(n) * normalOffsetScale;

        shadow_pos += vec4(normalOffset, 0.0);
    }

    shadow_coord = depth_vp_offset * model * shadow_pos;

    vec4 vpoint_mv = view * model * vec4(vpoint, 1.0);
    vpoint_MV = vpoint_mv;

    gl_Position = projection * vpoint_mv;
    distToCamera = gl_Position.z;
    uv = vtexcoord;
}
