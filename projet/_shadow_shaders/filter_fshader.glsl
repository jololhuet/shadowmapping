#version 330 core

in vec3 n;  // incoming normal;

out vec3 color;

uniform bool use_color;  // Use predefined color or texture?
uniform float texRatio;

uniform int light_type;

uniform float pcfSpread;
uniform float nIter = 16.0;


// For texturing
uniform sampler2D tex;
in vec2 uv;

// For coloring
uniform vec3 mesh_color;

uniform vec3 light_d;
uniform vec3 light_pos;

uniform float bias;
uniform bool usePolygonOffset;

uniform sampler2D shadow_map;
in vec4 shadow_coord;
in vec4 vpoint_MV;

float M_PI = 3.14159265358979323846;


// Poisson disk sample locations.
vec2 poisson_disk[16] = vec2[](
   vec2(-0.94201624, -0.39906216),
   vec2(0.94558609, -0.76890725),
   vec2(-0.094184101, -0.92938870),
   vec2(0.34495938, 0.29387760),
   vec2(-0.91588581, 0.45771432),
   vec2(-0.81544232, -0.87912464),
   vec2(-0.38277543, 0.27676845),
   vec2(0.97484398, 0.75648379),
   vec2(0.44323325, -0.97511554),
   vec2(0.53742981, -0.47373420),
   vec2(-0.26496911, -0.41893023),
   vec2(0.79197514, 0.19090188),
   vec2(-0.24188840, 0.99706507),
   vec2(-0.81409955, 0.91437590),
   vec2(0.19984126, 0.78641367),
   vec2(0.14383161, -0.14100790)
);

/*
 * Inspired by the perlin noise pseudo random
 */
float randNumber (vec4 pos, int i) {
    float bigger = 43758.5453;
    vec4 choosen = vec4(12.9898, 78.233, 45.16, 94.673);

    return fract(sin(dot(pos * (i+1) ,vec4(choosen))) * bigger);
}

float randomAngle (vec4 pos, int i) {
    return randNumber(pos, i) * 6.283285;
}

void main() {
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
    float shade = ambient_light + max(dot(normalize(n), normalize(light_dir)), 0.0);

    float shadow = 1.0;


    float spread = clamp(pcfSpread, 30, 5030);

    float randAngle = randomAngle(vpoint_MV, 16);
    float s = sin(randAngle);
    float c = cos(randAngle);

    for (int i = 0; i < int(nIter); i++) {
        int index = i;
//        int index = int(16.0*randNumber(shadow_coord, i))%16; // Will change with the camera moving
        vec2 rotatedOffset = vec2(poisson_disk[index].x * s, poisson_disk[index].y * c);

        if (texture2D(shadow_map, (shadow_coord.xy/shadow_coord.w + rotatedOffset/spread)).r < (shadow_coord.z - usedBias)/shadow_coord.w) {
            shadow -= (1.0/nIter);
        }
    }

    if (use_color) {
        color = shadow * shade * mesh_color;
    } else {
        color = shadow * shade * texture2D(tex, uv*texRatio).rgb;
    }
}
