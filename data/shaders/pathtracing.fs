
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
uniform float u_Samples;

vec2 randState = vec2(0.0);
float seed = 0.0;
out vec4 FragColor;

#define M_PI 3.141592653589

float InterleavedGradientNoise(vec2 xy) { return fract(52.9829189 * fract(xy.x * 0.06711056 + xy.y * 0.00583715)); }

float rand()
{
    float result = fract(sin((seed+1.0/u_Samples*InterleavedGradientNoise(gl_FragCoord.xy+1.0/u_Samples)) / 100.0 * dot(gl_FragCoord.xy, vec2(12.9898, 78.233))) * 43758.5453);
    seed += 1.0;
    return result;
}
float nextFloat(float a, float b) { 
	/*float rndx = InterleavedGradientNoise(gl_FragCoord.xy/textureSize(u_Texture, 0) + randState + vec2(u_Samples)); 
	randState += vec2(rndx, rndx); 
	float rndy = InterleavedGradientNoise(gl_FragCoord.xy/textureSize(u_Texture, 0) + randState + vec2(u_Samples)); 
	randState += vec2(rndy, rndy); 
	return mod(rndx - rndy, b); 
	*/
	return clamp((b-a)*rand(), a, b);
}
vec3 randomUnitVector()
{
	float z = nextFloat(0.0, 1.0) * 2.0 - 1.0;
	float a = nextFloat(0.0, 1.0) * 2.0 * M_PI;
	float r = sqrt(1.0 - z*z);
	float x = cos(a) * r;
	float y = sin(a) * r;
	return vec3(x, y, z);
}

const float uvFactor = (1.001 / 10.0);

void main()
{
	const float sunSolidAngle = 1e-5*1000.0 *19000.0*0.01;

    // generate primary ray
    vec2 ij = gl_FragCoord.xy / textureSize(u_Texture, 0).xy;
    vec3 ro = u_CameraPosition;
    vec3 rd = normalize(u_CameraLowerLeftCorner + ij.x*u_CameraHorizontal + ij.y*u_CameraVertical - ro);
	float mint = 0.001;
    float maxt = 10000.0;
    
    vec3 accum = vec3(0.0);
    vec3 emission = vec3(1.0);
    
    for (int level = 0; level < u_Depth; ++level)
    {
        float recT;
        vec3 recPoint;
        vec3 recNormal;
		float hitIndex = -1.0;

		// intersect world
        for (int obj = 0; obj < 10; ++obj)
        {
			float id = float(obj);
            vec4 data = texture(u_SpheresTexture, vec2(uvFactor*id, 0.0));
            vec3 center = data.xyz;
            float radius = data.w;
			float invRadius = 1.0 / data.w;
            vec3 to = ro - center;
            
			float B = dot(to, rd);
            float C = dot(to, to) - radius*radius;
            float discriminant = B*B - C;
            
            if (discriminant > 0.0) {
                float rootDiscriminant = sqrt(discriminant);
                float temp = (-B - rootDiscriminant);
                float temp2 = (-B + rootDiscriminant);
                if (temp < maxt && temp > mint) {
                    recPoint = ro + temp*rd;
                    recNormal = (recPoint - center) * invRadius; recT = temp;
					hitIndex = id;
					maxt = recT;
                }
				// this is required for transparent dielectric
                else if (temp2 < maxt && temp2 > mint) {
                    recPoint = ro + temp2*rd;
                    recNormal = (recPoint - center) * invRadius; recT = temp2;
					hitIndex = id;
					maxt = recT;
                }
            }
        }
        
		if (hitIndex >= 0) {
			//maxt = recT;
			vec4 material = texture(u_MaterialsTexture, vec2(uvFactor*hitIndex, 0.0));
			vec3 N = recNormal;
			ro = recPoint + N*0.001;
			mint = 0.001;
			maxt = 1.0;
			vec3 attenuation = material.bgr;
			if (material.w < 0.5) 
			{
				//float NdotL = -N.z;
				//attenuation *= max(NdotL, 0.0);
				
				// Diffuse BRDF (albedo/PI) with uniform Monte-Carlo 
				emission *= attenuation * 2.f;
				
				// Diffuse scattering
				rd = normalize(N + randomUnitVector());

				//ro = vec3(-1000.0, -1000.0, -1000.0);
				//rd = normalize(vec3(-1000.0, -1000.0, -1000.0));

				/*float sunLight = NdotL;
				if (sunLight > 0.0) {				
					vec3 direct = material.rgb / M_PI * sunLight * sunSolidAngle;
					accum += emission * direct;
				}*/
				
			}
			/*
			else if (material.w < 1.5) {
				emission *= attenuation;
				ro = recPoint;
				rd = reflect(rd, N);
				//if (dot(recNormal, rd) < 0.0) {
				//	ro = vec3(-10000.0, -10000.0, -10000.0);
				//break;
				//}
			}
			else {
				ro = recPoint;
				rd = refract(rd, N, 1.55);
				attenuation = vec3(1.0);
			}
			*/

		}
        else {
            float t = 0.5 * (rd.y + 1.0);
            vec3 skyLight = vec3(1.0 - t) + vec3(0.5, 0.7, 1.0) * (t);
            accum += emission * skyLight;
            break;
        }
    }
    
    vec3 source = vec3(u_Samples);// + texture(u_Texture, ij).rgb;
    FragColor = vec4(source, 1.0);
}
