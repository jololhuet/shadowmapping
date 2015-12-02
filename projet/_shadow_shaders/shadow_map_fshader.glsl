#version 330 core
uniform bool isVsm = false;

uniform float bias = -1.0;

in vec4 vpoint_mv;

out vec4 color;

void main()
{
    if (isVsm) {
        float depth = vpoint_mv.z / vpoint_mv.w;
        depth = depth * 0.5 + 0.5;

        float moment1 = depth;
        float moment2 = depth * depth;

        // Adjusting moments (this is sort of bias per pixel) using partial derivative
        if (bias == -1){
            float dx = dFdx(depth);
            float dy = dFdy(depth);
            moment2 += 0.25*(dx*dx+dy*dy);
        } else {
            float dx = bias;
            float dy = bias;
            moment2 += 0.25*(dx*dx+dy*dy);
        }

        color = vec4(moment1, moment2, 0.0, 0.0);
    }
}
