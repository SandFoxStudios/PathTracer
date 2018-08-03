
precision mediump float;

in vec2 v_TexCoords;

uniform sampler2D u_Texture;

uniform float u_InvNumSamples;

out vec4 FragColor;

void main()
{
    vec3 source = /*vec3(u_InvNumSamples) **/ texture(u_Texture, v_TexCoords).rgb;
    source = pow(source, vec3(0.454545));
    FragColor = vec4(source, 1.0);
}
