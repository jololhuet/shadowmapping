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

uniform sampler2D shadow_map;
in vec4 shadow_coord;
in vec4 vpoint_MV;


void main() {
    vec3 light_dir;
    if (light_type == 0) {
        light_dir = light_d;
    } else {
        light_dir = normalize(light_pos-vpoint_MV.xyz);
    }

    float ambient_light = 0.1;
    float shade = ambient_light + max(dot(normalize(n), normalize(light_dir)), 0.0);

    float shadow = 1.0;

    if (texture2D(shadow_map, (shadow_coord.xy / shadow_coord.w)).r < shadow_coord.z / shadow_coord.w - 0.001) {
        shadow = 0.2;
    }

    if (use_color) {
        color = shadow * shade * mesh_color;
    } else {
        color = shadow * shade * texture2D(tex, uv*texRatio).rgb;
    }
}
