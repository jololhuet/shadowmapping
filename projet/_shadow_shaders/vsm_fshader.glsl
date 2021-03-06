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

uniform bool useCsm;

uniform sampler2D vsm_shadow_map0;
uniform sampler2D vsm_shadow_map1;
uniform sampler2D vsm_shadow_map2;
uniform sampler2D vsm_shadow_map3;
uniform sampler2D vsm_shadow_map4;
uniform sampler2D vsm_shadow_map5;
uniform sampler2D vsm_shadow_map6;
uniform sampler2D vsm_shadow_map7;

uniform float farPlane[8];

in vec4 shadow_coord;
in vec4 vpoint_MV;
in float distToCamera;

// Uniforms for blur
uniform float blurRadius = 0.01;
uniform int blurSteps = 4;

// http://www.ogre3d.org/forums/viewtopic.php?f=5&t=50786
vec4 blurTex2D(sampler2D map, vec2 uv_coord, float radius, int steps)
   {
      float stepSize = 2.0 * radius / (1.0*steps);
      uv_coord -= vec2(radius);

      vec4 total = vec4(0, 0, 0, 0);

      for (int x = 0; x < steps; ++x)
         for (int y = 0; y < steps; ++y)
            total += texture2D(map, vec2(uv_coord + vec2(x * stepSize, y * stepSize)));

      return total / (steps * steps);
   }

float chebyshevUpperBound(vec4 coords, float minVar, sampler2D vsm_shadow_map)
{
    // We retrive the two moments previously stored (depth and depth*depth)
    vec2 moments = blurTex2D(vsm_shadow_map, coords.xy/coords.w, blurRadius, blurSteps).xy;

    // Surface is fully lit. as the current fragment is before the light occluder
    if (coords.z/coords.w <= moments.x)
            return 1.0 ;

    // The fragment is either in shadow or penumbra. We now use chebyshev's upperBound to check
    // How likely this pixel is to be lit (p_max)
    float variance = min(max(moments.y - (moments.x*moments.x), minVar), 1.0);
//    variance = max(variance, minVar);

    float distance = coords.z / coords.w;

    float d = distance - moments.x;
    float p_max = variance / (variance + (d*d));

    float amount = 0.5f;
    p_max =  clamp((p_max - amount) / (1.0f - amount), 0.01f, 1.0f); // standard light bleeding fix

    return p_max;
}

vec3 proceed (sampler2D vsm_shadow_map) {
    vec3 light_dir;
    if (light_type == 0) {
        light_dir = light_d;
    } else {
        light_dir = normalize(light_pos-vpoint_MV.xyz);
    }

    float ambient_light = 0.1;
    float shade = ambient_light + max(dot(normalize(n), normalize(light_dir)), 0.0);

    float shadow = chebyshevUpperBound(shadow_coord, 0.0001f, vsm_shadow_map);

    if (use_color) {
        return shadow * shade * mesh_color;
    } else {
        return shadow * shade * texture2D(tex, uv*texRatio).rgb;
    }
}

void main() {
    if (useCsm) {
        if (distToCamera < farPlane[0]) {
            color = proceed(vsm_shadow_map0);
        } else if (distToCamera < farPlane[1]) {
            color = proceed(vsm_shadow_map1);
        } else if (distToCamera < farPlane[2]) {
            color = proceed(vsm_shadow_map2);
        } else if (distToCamera < farPlane[3]) {
            color = proceed(vsm_shadow_map3);
        } else if (distToCamera < farPlane[4]) {
            color = proceed(vsm_shadow_map4);
        } else if (distToCamera < farPlane[5]) {
            color = proceed(vsm_shadow_map5);
        } else if (distToCamera < farPlane[6]) {
            color = proceed(vsm_shadow_map6);
        } else if (distToCamera < farPlane[7]) {
            color = proceed(vsm_shadow_map7);
        } else {
            vec3(0);
        }
    } else {
        color = proceed(vsm_shadow_map0);
    }
}
