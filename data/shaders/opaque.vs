
#extension GL_ARB_explicit_attrib_location : enable

layout(location=0) in vec3 a_Position;
layout(location=1) in vec3 a_Normal;
layout(location=2) in vec3 a_TexCoords;
layout(location=3) in vec4 a_Color;

uniform mat4 u_ViewMatrix;
uniform mat4 u_ProjectionMatrix;

out vec3 v_Position;
out vec3 v_Normal;
out vec3 v_TexCoords;
out vec4 v_Color;

// Render a screen-space quad
void main(void)
{
    v_TexCoords = a_TexCoords;
    v_Normal = mat3(u_ViewMatrix) * a_Normal;
    v_Color = vec4(/*GammaToLinearSpace*/(a_Color.rgb), a_Color.a);
    
    gl_Position = u_ProjectionMatrix * u_ViewMatrix * vec4(a_Position, 1.0);
}
