
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

#define M_PI 3.141592653589

float InterleavedGradientNoise(vec2 xy) { return fract(52.9829189 * fract(xy.x * 0.06711056 + xy.y * 0.00583715)); }

/*vec3 randomUnitVector()
{
	float z = nextFloat(0.0, 1.0) * 2.0 - 1.0;
	float a = nextFloat(0.0, 1.0) * 2.0 * M_PI;
	float r = sqrtf(1.0 - z*z);
	float x = cos(a) * r;
	float y = sin(a) * r;
	return vec3(x, y, z);
}*/

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
		float hitIndex = -1.0;
        bool hit = false;
        
		const float gradient = 1.0 / 10.0;

		// intersect world
        for (int obj = 0; obj < 10; ++obj)
        {
			float id = float(obj);
            vec4 data = texture(u_SpheresTexture, vec2(gradient*(id), 0.0));
            vec3 center = data.xyz;
            float radius = data.w;
            vec3 to = ro - center;
            
            float B = dot(to, rd);
            float C = dot(to, to) - radius*radius;
            float discriminant = B*B - A*C;
            
            if (discriminant > 0.0) {
                float rootDiscriminant = sqrt(discriminant);
                float temp = (-B - rootDiscriminant) * invA;
                float temp2 = (-B + rootDiscriminant) * invA;
                if (temp < maxt && temp > mint) {
                    recPoint = ro + temp*rd;
                    recNormal = (recPoint - center) * (1.0/radius); recT = temp;
                    hit = true;
					hitIndex = id;
                }
                else if (temp2 < maxt && temp2 > mint) {
                    recPoint = ro + temp2*rd;
                    recNormal = (recPoint - center) * (1.0/radius); recT = temp2;
                    hit = true;
					hitIndex = id;
                }
            }
        }
        
		if (hitIndex >= 0.0) {
			vec3 attenuation = vec3(1.0);
			vec4 material = texture(u_MaterialsTexture, vec2(gradient*(hitIndex), 0.0));
			attenuation = material.rgb / M_PI;// * max(-normalize(recNormal).z, 0.0);
			emission *= attenuation;

			ro = recPoint;
			rd = recNormal;
			mint = 0.001;
			maxt = 100.0;

			if (material.w < 0.5 && -recNormal.z > 0)
				accum += emission;
		}
        else {
            float t = 0.5 * (rd.y + 1.0);
            vec3 skyLight = vec3(1.0)*(1.0 - t) + vec3(0.5, 0.7, 1.0) * (t);
            accum += emission * skyLight;
            break;
        }
    }
    
    vec3 source = accum;// + texture(u_Texture, v_TexCoords).rgb;
    FragColor = vec4(source, 1.0);
}
