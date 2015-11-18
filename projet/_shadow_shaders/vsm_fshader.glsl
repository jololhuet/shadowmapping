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

uniform sampler2D vsm_shadow_map;

in vec4 shadow_coord;
in vec4 vpoint_MV;

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

float chebyshevUpperBound(vec4 coords, float minVar)
{
    // We retrive the two moments previously stored (depth and depth*depth)
    vec2 moments = blurTex2D(vsm_shadow_map,coords.xy/coords.w, blurRadius, blurSteps).xy;

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


void main() {
    vec3 light_dir;
    if (light_type == 0) {
        light_dir = light_d;
    } else {
        light_dir = normalize(light_pos-vpoint_MV.xyz);
    }

    float ambient_light = 0.1;
    float shade = ambient_light + max(dot(normalize(n), normalize(light_dir)), 0.0);

    float shadow = chebyshevUpperBound(shadow_coord, 0.0001f);

    if (use_color) {
        color = shadow * shade * mesh_color;
    } else {
        color = shadow * shade * texture2D(tex, uv*texRatio).rgb;
    }
}
