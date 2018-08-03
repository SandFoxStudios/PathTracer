
const vec2 quad[4] = vec2[]
(
 vec2(-1.0,  1.0),
 vec2(-1.0, -1.0),
 vec2( 1.0,  1.0),
 vec2( 1.0, -1.0)
);

out vec3 v_Position;
out vec2 v_TexCoords;

// Render a screen-space quad
void main(void)
{
    v_Position = vec3(quad[gl_VertexID], -1.0);
    v_TexCoords = v_Position.xy * 0.5 + 0.5;
    
    vec3 a_Position = v_Position;
    gl_Position = vec4(a_Position, 1.0);
}
