
precision highp float;

in vec2 v_TexCoords;

uniform vec3 u_CameraPosition;
uniform vec3 u_CameraLowerLeftCorner;
uniform vec3 u_CameraHorizontal;
uniform vec3 u_CameraVertical;

uniform sampler2D u_Texture;
uniform sampler2D u_SpheresTexture;
uniform sampler2D u_MaterialsTexture;

uniform int u_Depth;

out vec4 FragColor;

void main()
{
    // generate primary ray
    vec2 ij = gl_FragCoord.xy / textureSize(u_Texture, 0).xy;
    vec3 ro = u_CameraPosition;
    vec3 rd = normalize(u_CameraLowerLeftCorner + ij.x*u_CameraHorizontal + ij.y*u_CameraVertical - ro);

    
    vec3 accum = vec3(0.0);
    vec3 emission = vec3(1.0);
    
    float A = dot(rd, rd);
    float invA = 1.0 / A;
    
    for (int level = 0; level < u_Depth; ++level)
    {
        float mint = 0.001;
        float maxt = 100.0;
        float recT;
        vec3 recPoint;
        vec3 recNormal;
        bool hit = false;
        
        for (int obj = 0; obj < 10; ++obj)
        {
            vec4 data = texture(u_SpheresTexture, vec2(1.0/(float(obj)+0.0), 0.0));
            vec3 center = data.xyz;
            float radius = data.w;
            vec3 to = ro - center;
            
            float B = dot(to, rd);
            float C = dot(to, to) - radius*radius;
            float discriminant = B*B - A*C;
            
            if (discriminant > 0.0) {
                float rootDiscriminant = sqrt(discriminant);
                float temp = (-B - rootDiscriminant) * invA;
                if (temp < maxt && temp > mint) {
                    recPoint = ro + temp*rd;
                    recNormal = (recPoint - center) * (1.0/radius); recT = temp;
                    hit = true;
                }
                float temp2 = (-B + rootDiscriminant) * invA;
                if (temp2 < maxt && temp2 > mint) {
                    recPoint = ro + temp2*rd;
                    recNormal = (recPoint - center) * (1.0/radius); recT = temp2;
                    hit = true;
                }
            }
        }
        
        if (!hit) {
            float t = 0.5 * (rd.y + 1.0);
            vec3 skyLight = vec3(1.0)*(1.0 - t) + vec3(0.5, 0.7, 1.0) * (t);
            accum = skyLight;
            break;
        }
    }
    
    vec3 source = accum;// + texture(u_Texture, v_TexCoords).rgb;
    FragColor = vec4(source, 1.0);
}
