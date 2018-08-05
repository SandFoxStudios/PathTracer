
precision highp float;

in vec2 v_TexCoords;

uniform sampler2D u_Texture;

uniform float u_InvNumSamples;

out vec4 FragColor;

void main()
{
	vec2 ij = gl_FragCoord.xy / textureSize(u_Texture, 0).xy;
    vec3 source = texture(u_Texture, ij).rgb / vec3(u_InvNumSamples);
	//source /= (source + vec3(1.0));
	source = pow(source, vec3(0.454545));
    FragColor = vec4(source, 1.0);
}
