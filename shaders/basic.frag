#version 330 core
out vec4 FragColor;

in float Color;
in vec2 Pos;

uniform float u_radius;

vec3 colormap(float t) {
    t = clamp(t, 0.0, 1.0);
    return vec3(
        clamp(1.5 - abs(4.0 * t - 3.0), 0.0, 1.0), // Red
        clamp(1.5 - abs(4.0 * t - 2.0), 0.0, 1.0), // Green
        clamp(1.5 - abs(4.0 * t - 1.0), 0.0, 1.0)  // Blue
    );
}

void main()
{
    float dist = length(Pos);
    vec3 rgb = colormap(Color);

    if (dist > u_radius) {
        FragColor = vec4(rgb, 1.0);
    } else {
        FragColor = vec4(1.0); // white
    }
}