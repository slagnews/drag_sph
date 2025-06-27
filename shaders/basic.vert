#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in float aColor;

out float Color;
out vec2 Pos;

void main()
{
    gl_Position = vec4(aPos, 0.0, 1.0);
    gl_PointSize = 5.0;
    Pos = aPos;
    Color = aColor;
}