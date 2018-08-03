
precision mediump float;

in vec3 v_Normal;
in vec3 v_TexCoords;
in vec4 v_Color;

uniform sampler2D u_Texture;
uniform vec4 u_MaterialColor;

out vec4 FragColor;

void main()
{
    vec4 texColor = textureProj(u_Texture, v_TexCoords);

    float alpha = texColor.a * v_Color.a * u_MaterialColor.a;

    vec3 color = texColor.rgb * v_Color.rgb * u_MaterialColor.rgb;
    FragColor = vec4(color.rgb * alpha, alpha);
}
