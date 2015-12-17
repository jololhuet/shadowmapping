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

uniform bool showFrustSplit;

in vec4 shadow_coord;
in vec4 vpoint_MV;
in float distToCamera;

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
    if (distToCamera < farPlane[0]) {
        if (!showFrustSplit){
            color = proceed(shadow_map0);
        } else {
            color = vec3(1, 0, 0);
        }
    } else if (distToCamera < farPlane[1]) {
        if (!showFrustSplit){
            color = proceed(shadow_map1);
        } else {
            color = vec3(0, 1, 0);
        }
    } else if (distToCamera < farPlane[2]) {
        if (!showFrustSplit){
            color = proceed(shadow_map2);
        } else {
            color = vec3(0, 0, 1);
        }
    } else if (distToCamera < farPlane[3]) {
        if (!showFrustSplit){
            color = proceed(shadow_map3);
        } else {
            color = vec3(1, 1, 0);
        }
    } else if (distToCamera < farPlane[4]) {
        if (!showFrustSplit){
            color = proceed(shadow_map4);
        } else {
            color = vec3(1, 0, 1);
        }
    } else if (distToCamera < farPlane[5]) {
        if (!showFrustSplit){
            color = proceed(shadow_map5);
        } else {
            color = vec3(0, 1, 1);
        }
    } else if (distToCamera < farPlane[6]) {
        if (!showFrustSplit){
            color = proceed(shadow_map6);
        } else {
            color = vec3(1, 0, 1);
        }
    } else if (distToCamera < farPlane[7]) {
        if (!showFrustSplit){
            color = proceed(shadow_map7);
        } else {
            color = vec3(1, 1, 0);
        }
    } else {
        vec3(0);
    }
}
