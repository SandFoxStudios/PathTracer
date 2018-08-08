
//precision highp float;

in vec2 v_TexCoords;

uniform vec3 u_CameraPosition;
uniform vec3 u_CameraLowerLeftCorner;
uniform vec3 u_CameraHorizontal;
uniform vec3 u_CameraVertical;

uniform sampler2D u_Texture;
uniform sampler2D u_SpheresTexture;
uniform sampler2DRect u_VerticesTexture;
uniform usampler2DRect u_IndicesTexture;
uniform sampler2D u_MaterialsTexture;

uniform int u_NumTriangles;
uniform int u_Depth;
uniform float u_Samples;

vec2 randState = vec2(0.0);
float seed = 0.0;
out vec4 FragColor;

#define M_PI 3.141592653589

float InterleavedGradientNoise(vec2 xy) { return fract(52.9829189 * fract(xy.x * 0.06711056 + xy.y * 0.00583715)); }

float rand()
{
	float result = InterleavedGradientNoise(gl_FragCoord.xy + randState.x + u_Samples);
	randState += vec2(result, result*2.0);
    //float result = fract(sin(seed / 100.0 * dot(gl_FragCoord.xy, vec2(12.9898, 78.233))) * 43758.5453);
    //seed += 1.0;//result;
    return result;
}
float rand01() { 
	/*float rndx = InterleavedGradientNoise(gl_FragCoord.xy/textureSize(u_Texture, 0) + randState + vec2(u_Samples)); 
	randState += vec2(rndx, rndx); 
	*/
	return mod(rand(), 1.0);
}

mat3 GetTangentSpace(vec3 normal)
{
    // Choose a helper vector for the cross product
    vec3 helper = vec3(1.0, 0.0, 0.0);
    if (abs(normal.x) > 0.99)
        helper = vec3(0.0, 0.0, 1.0);

    // Generate vectors
    vec3 tangent = normalize(cross(normal, helper));
    vec3 binormal = normalize(cross(normal, tangent));
    return mat3(tangent, normal, binormal);
}

vec3 randomUnitVector()
{
	float z = rand01() * 2.0 - 1.0;
	float r = sqrt(1.0 - z*z);
	float a = rand01() * 2.0 * M_PI;
	float x = cos(a) * r;
	float y = sin(a) * r;
	return vec3(x, y, z);
}

vec3 randomHemisphere(float r1, float distrib)
{
    float cosTheta = pow(r1, 1.0 / (distrib+1.0));
	float sinTheta = sqrt(max(0.0, 1.0 - cosTheta*cosTheta));
	float phi = rand01() * 2.0 * M_PI;
	return vec3(cos(phi) * sinTheta,
				cosTheta,
                sin(phi) * sinTheta
				);
}

const float uvFactor = (1.001 / 10.0);
const float sunSolidAngle = 1e-5*1000.0 *19000.0*0.01;
const float maxt = 1.0e38;
const float mint = 0.00;

float IntersectWorld(vec3 ro, vec3 rd, inout vec3 recPoint, inout vec3 recNormal, inout float recT)
{
	float hitIndex = -1.0;
    
	for (int obj = 0; obj < 10; ++obj)
    {
		float id = float(obj);
        vec4 data = texture(u_SpheresTexture, vec2(uvFactor*id, 0.0));
        vec3 center = data.xyz;
        float radius = data.w;
		float invRadius = 1.0 / data.w;
        vec3 to = ro - center;
            
		float B = -dot(to, rd);
        float C = dot(to, to) - radius*radius;
        float discriminant = B*B - C;
            
        if (discriminant > 0.0) {
            float rootDiscriminant = sqrt(discriminant);
            float temp = (B - rootDiscriminant);
			temp = temp > 0.0 ? temp : (B + rootDiscriminant);
            //float temp2 = (B + rootDiscriminant);
            if (temp > mint && temp < recT) {
                recPoint = ro + temp*rd;
                recNormal = (recPoint - center)*invRadius; 
				recT = temp;
				hitIndex = id;
            }
			// this is required for transparent dielectric
            /*else if (temp2 > mint && temp2 < recT) {
                recPoint = ro + temp2*rd;
                recNormal = normalize(recPoint - center); recT = temp2;
				hitIndex = id;
            }*/
        }
    }
    
	return hitIndex;
}

const float texFactor = 1.0;// / 1024.0;

float IntersectMesh(vec3 ro, vec3 rd, inout vec3 recPoint, inout vec3 recNormal, inout float recT)
{
    float hitIndex = -1.0;
    
    for (int i = 0; i < u_NumTriangles; i++)
    {
        uvec3 indices = texelFetch(u_IndicesTexture, ivec2(i, 0)).stp;
        uvec3 vind = uvec3(indices*2);
        vec3 a = texelFetch(u_VerticesTexture, ivec2(vind.s, 0)).xyz;
        vec3 b = texelFetch(u_VerticesTexture, ivec2(vind.t, 0)).xyz;
        vec3 c = texelFetch(u_VerticesTexture, ivec2(vind.p, 0)).xyz;
        vec3 e1 = b - a;
        vec3 e2 = c - a;

        vec3 q = cross(rd, e2);
        float det = dot(e1, q);
        
        if (det >= 0.0001)
        {
            vec3 ap = (ro - a);
            vec3 r = cross(ap, e1);
            
            float u = dot(ap, q);
            float v = dot(rd, r);
            
            
            if (u >= 0.0 && u <= det && v >= 0.0 && (u + v) <= det)
            {
                vec3 uv;
                uvec3 nind = vind+ivec3(1);
                vec3 na = texelFetch(u_VerticesTexture, ivec2(nind.s, 0)).xyz;
                vec3 nb = texelFetch(u_VerticesTexture, ivec2(nind.t, 0)).xyz;
                vec3 nc = texelFetch(u_VerticesTexture, ivec2(nind.p, 0)).xyz;
                
                float ood = 1.0 / det;
                float temp = dot(e2, r) * ood;
                if (temp > mint && temp < recT) {
                    recT = temp;
                    recPoint = ro + temp*rd;
                    uv.x = u*ood; uv.y = v*ood; uv.z = 1.f - uv.x - uv.y;
                    recNormal = normalize(na * uv.z + nb * uv.x + nc * uv.y);
                    //recNormal = normalize(cross(e1, e2));
                    hitIndex = 0.0;
                }
            }
        }
    }
    
    return hitIndex;
}

float FurnaceHitTest(vec3 ro, vec3 rd, inout vec3 recPoint, inout vec3 recNormal, inout float recT)
{
    float hitIndex = -1.0;
    
    vec3 center = vec3(0.0, 0.0, 1.41);
    float radius = 0.5;
    float invRadius = 1.0 / radius;
    vec3 to = ro - center;
    
    float B = -dot(to, rd);
    float C = dot(to, to) - radius*radius;
    float discriminant = B*B - C;
    
    if (discriminant > 0.0) {
        float rootDiscriminant = sqrt(discriminant);
        float temp = (B - rootDiscriminant);
        temp = temp > 0.0 ? temp : (B + rootDiscriminant);
        //float temp2 = (B + rootDiscriminant);
        if (temp > mint && temp < recT) {
            recPoint = ro + temp*rd;
            recNormal = (recPoint - center)*invRadius;
            recT = temp;
            hitIndex = 100;
        }
        // this is required for transparent dielectric
        /*else if (temp2 > mint && temp2 < recT) {
         recPoint = ro + temp2*rd;
         recNormal = normalize(recPoint - center); recT = temp2;
         hitIndex = id;
         }*/
    }
    
    return hitIndex;
}

float saturate(float x) { return clamp(x, 0.0, 1.0); }
// direct lighting Fresnel/Schlick (technically non metals ?)
float SphericalGaussianApprox(float CosX, float ModifiedSpecularPower)
{
    return exp2(ModifiedSpecularPower * CosX - ModifiedSpecularPower);
}
#define OneOnLN2_x6 8.656170 // == 1/ln(2) * 6   (6 is SpecularPower of 5 + 1)
vec3 FresnelSchlick(vec3 SpecularColor, vec3 E, vec3 H)
{
    // In this case SphericalGaussianApprox(1.0f - saturate(dot(E, H)), OneOnLN2_x6) is equal to exp2(-OneOnLN2_x6 * x)
    return SpecularColor + (vec3(1.0) - SpecularColor) * exp2(-OneOnLN2_x6 * saturate(dot(E, H)));
}

vec3 Shade(inout vec3 attenuation, inout vec3 ro, inout vec3 rd, float hitIndex, vec3 P, vec3 N)
{
	if (hitIndex >= 0.0) 
	{
        vec4 material = texture(u_MaterialsTexture, vec2(uvFactor*hitIndex, 0.0));
        material = vec4(0.5);

        // white furnace Monte-Carlo test
        /*vec3 albedo = material.rgb / M_PI;
        vec3 indirect = vec3(0.0);
        float recT = maxt;
        int count = 256;
        float pdf = 1.0 / (2.0*M_PI);
        for (int n = 0; n < count; n++)
        {
            float r1 = rand01();
            rd = (GetTangentSpace(N)*randomHemisphere(r1));
            float p = r1/pdf;
            if (FurnaceTest(ro, rd, P, N, recT) >= 0)
                indirect += vec3(p, 0.0, 0.0);
            else {
                float t = 0.5 * (rd.y + 1.0);
                vec3 skyLight = vec3(1.0 - t) + vec3(0.5, 0.7, 1.0) * (t);
                indirect += vec3(p) * skyLight;
            }
        }
        indirect /= vec3(float(count));
        attenuation *= (indirect) * albedo;*/
		// test glossy plastic ---
		vec3 specularity = vec3(0.94);
        const float roughness = 0.0;
		const float smoothness = 1.0 - roughness;
        const float alpha = pow(1024.0, smoothness*smoothness);//15.0;
		vec3 albedo = material.bgr;
		
		vec3 R = reflect(rd, N);
		// perfect reflection
		//rd = R;
        //ro = P + rd*0.0001;
		//attenuation *= albedo;
        //return vec3(0.0);
        //return albedo * max(-N.z, 0.0);//dot(N, rd)
		
        float r1 = rand01();
        
		// Lambert ...
        //float distrib = 0.0; // 0 if uniform, 1 if cosine-weighted
        //rd = (GetTangentSpace(N)*randomHemisphere(r1, distrib));
        //rd = normalize(N + randomUnitVector());
		//rd = mix(R, rd, roughness);
        //ro = P + rd*0.0001;
		// ...uniform BRDF (pdf = 1.0 / 2PI)
        //vec3 lambert = albedo * 2.0 * clamp(dot(N, rd), 0.0, 1.0);
		// ...importance BRDF (pdf = cos@ / PI)
		//vec3 lambert = (albedo);
		//attenuation *= lambert;

        // Phong ...
        float distrib = alpha;//0.0; // 0 if uniform, alpha if cosine-weighted
        rd = (GetTangentSpace(R)*randomHemisphere(r1, distrib));
		// modified-Phong uniform BRDF (pdf = 1.0 / 2PI)
        //vec3 diffuse = vec3(0.0);//min(vec3(1.0) - specularity, albedo) * 2.0;
		//vec3 specular = specularity * (alpha + 2.0) * pow(max(dot(rd, R), 0.0), alpha);
		//attenuation *= (diffuse + specular) * clamp(dot(N, rd), 0.0, 1.0);
        // modified-Phong uniform BRDF (pdf = ((n+1)/2PI)cos@^n / PI)
        vec3 diffuse = min(vec3(1.0) - specularity, albedo);
        float f = (alpha + 2.0)/(alpha+1.0);
        vec3 specular = FresnelSchlick(specularity, N, rd) * specularity * clamp(dot(N, rd) * f, 0.0, 1.0);
        attenuation *= (diffuse + specular);
		
		// ... with perfect reflection (russian roulette)
		/*vec3 diffuse = min(vec3(1.0) - specularity, albedo);
		float specChance = dot(specularity, vec3(1.0/3.0));
		float diffChance = dot(diffuse, vec3(1.0/3.0));
		float sum = specChance + diffChance;
		specChance /= sum;
		diffChance /= sum;
		float roulette = mod(rand(), 1.0);
		if (roulette < specChance) 
		{
			rd = reflect(rd, N);
			float factor = (alpha + 2.0) / (alpha + 1.0);
			float NdotR = clamp(dot(N, rd) * factor, 0.0, 1.0);
				
			attenuation *= (1.0/specChance) * specularity * NdotR;
				
		}
		else {
			//float NdotR = clamp(dot(N, rd), 0.0, 1.0);
			attenuation *= (1.0/diffChance) * diffuse;// * 2.0 * NdotR;
		}*/
		// ---

		/*if (material.w < 0.5) 
		{
			// Diffuse scattering
			rd = normalize(N + randomUnitVector());

			float NdotL =  max(-N.z, 0.0);
			//attenuation *= max(NdotL, 0.0);
				
			// Diffuse BRDF (albedo/PI) with uniform Monte-Carlo 
			attenuation = material.bgr * 2.f;

			//ro = vec3(-1000.0, -1000.0, -1000.0);
			//rd = normalize(vec3(-1000.0, -1000.0, -1000.0));

			//float sunLight = NdotL;
			/if (sunLight > 0.0) {				
			//	vec3 direct = material.rgb / M_PI * sunLight * sunSolidAngle;
			//	accum += energy * direct;
			//}

			// outgoing
			emission *= albedo * -rd.z;
		}
			
		else if (material.w < 1.5) {
			//attenuation = material.bgr;
				
			attenuation = vec3(0.0);
			// outgoing
			emission *= material.bgr;
				
			rd = normalize(reflect(rd, N) + 0.01*randomUnitVector());
			//if (dot(recNormal, rd) < 0.0) {
			//	ro = vec3(-10000.0, -10000.0, -10000.0);
			//break;
			//}
		}
		else {
			rd = refract(rd, N, 1.55);
			attenuation = vec3(0.0);
			emission = vec3(0.0);
		}*/
			
		return vec3(r1);
	}
    else 
	{
        float t = 0.5 * (rd.y + 1.0);
        vec3 skyLight = vec3(1.0 - t) + vec3(0.5, 0.7, 1.0) * (t);
		//skyLight *= attenuation;
		//attenuation = vec3(0.0);
		
		// outgoing
        return vec3(1.0);//skyLight;
    }
}



void main()
{
	vec3 accum = vec3(0.0);
	
    // generate primary ray
    vec2 ij = (gl_FragCoord.xy /*+ (rand01()-0.5)*/) / textureSize(u_Texture, 0).xy;
    vec3 ro = u_CameraPosition;
    vec3 rd = normalize(u_CameraLowerLeftCorner + ij.x*u_CameraHorizontal + ij.y*u_CameraVertical - ro);
    
    
    vec3 energy = vec3(1.0);
    vec3 prob = vec3(1.0);
    for (int level = 0; level < u_Depth; ++level)
    {
		float recT = maxt;
        vec3 recPoint;
        vec3 recNormal;

        //float hitIndex = FurnaceHitTest(ro, rd, recPoint, recNormal, recT);
		float hitIndex = IntersectWorld(ro, rd, recPoint, recNormal, recT);
        //float hitIndex = IntersectMesh(ro, rd, recPoint, recNormal, recT);
        
	    if (hitIndex>=0) {
            prob = Shade(energy, ro, rd, hitIndex, recPoint, recNormal);
		}
				
		if (!any(bvec3(energy)) || hitIndex<0.0) {
			break;
		}
    }

	// ambient
	float t = 0.5 * (rd.y + 1.0);
    vec3 skyLight = vec3(1.0 - t) + vec3(0.5, 0.7, 1.0) * (t);
    accum += energy * /*prob **/ skyLight;

    //vec3 source = vec3(accum) + texture(u_Texture, ij).rgb;
	// incremental averaging
	vec3 source = mix(texture(u_Texture, ij).rgb, accum, 1.0/u_Samples);
    FragColor = vec4(source, 1.0);
}
