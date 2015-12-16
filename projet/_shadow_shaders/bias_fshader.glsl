#version 330 core

in vec3 n;  // incoming normal;

out vec3 color;

uniform sampler2D normalTex;
uniform bool use_color;  // Use predefined color or texture?
uniform float texRatio;
uniform mat4 model;

// For texturing
uniform sampler2D tex;
in vec2 uv;

uniform vec3 light_d;
uniform vec3 light_pos;

uniform int light_type;

uniform float shadowDarkness;

// For coloring
uniform vec3 mesh_color;

uniform float bias;
uniform bool usePolygonOffset;

uniform bool useCsm;

uniform sampler2D shadow_map0;
uniform sampler2D shadow_map1;
uniform sampler2D shadow_map2;
uniform sampler2D shadow_map3;
uniform sampler2D shadow_map4;
uniform sampler2D shadow_map5;
uniform sampler2D shadow_map6;
uniform sampler2D shadow_map7;

uniform float farPlane[8];

in vec4 shadow_coord;
in vec4 vpoint_MV;

vec3 proceed (sampler2D shadow_map) {
    vec3 light_dir;
    if (light_type == 0) {
        light_dir = light_d;
    } else {
        light_dir = normalize(light_pos-vpoint_MV.xyz);
    }

    float usedBias = 0;

    if (!usePolygonOffset) {
        float cosTheta = normalize(dot(n, light_dir));
        usedBias = bias * tan(acos(cosTheta));
        usedBias = clamp(bias, 0.001, 20);
    }

    float ambient_light = 0.1;
    float shade = ambient_light + dot(normalize(n), (light_dir));

    float shadow = 1.0;  // shading factor from the shadow (1.0 = no shadow, 0.0 = all dark)

    if (texture2D(shadow_map0, (shadow_coord.xy / shadow_coord.w)).r < (shadow_coord.z-usedBias) / shadow_coord.w) {
        shadow = clamp(shadowDarkness, 0, 0.95);
    }

    if (use_color) {
        return shadow * shade * mesh_color;
    } else {
        return shadow * shade * texture2D(tex, uv*texRatio).rgb;
    }
}

void main() {
    if (vpoint_MV.z < farPlane[0]) {
        color = proceed(shadow_map0);
    } else if (vpoint_MV.z < farPlane[1]) {
        color = proceed(shadow_map1);
    } else if (vpoint_MV.z < farPlane[2]) {
        color = proceed(shadow_map2);
    } else if (vpoint_MV.z < farPlane[3]) {
        color = proceed(shadow_map3);
    } else if (vpoint_MV.z < farPlane[4]) {
        color = proceed(shadow_map4);
    } else if (vpoint_MV.z < farPlane[5]) {
        color = proceed(shadow_map5);
    } else if (vpoint_MV.z < farPlane[6]) {
        color = proceed(shadow_map6);
    } else if (vpoint_MV.z < farPlane[7]) {
        color = proceed(shadow_map7);
    } else {

    }
}
