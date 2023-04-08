#version 330

in highp vec3 v_vertex;
in highp vec3 v_normal;
in highp vec2 v_texcoord;
in highp float v_material;

out highp vec4 fragColor;

uniform highp vec3 lightPos;
uniform sampler2D textures[16];

void main() {
    highp vec3 L = normalize(lightPos - v_vertex);
    highp float NL = max(dot(normalize(v_normal), L), 0.0);

//   highp vec3 color = vec3(0.39, 1.0, 0.0);
//   highp vec3 col = clamp(color * 0.2 + color * 0.8 * NL, 0.0, 1.0);
//   fragColor = vec4(col, 1.0);

    int i;
    int material = int(v_material);
    highp vec4 mixTexture = vec4(0.0, 0.0, 0.0, 0.0);
    highp vec4 tex;
    for (i = 0; i < 16; ++i) {
      if ((material & (1 << i)) != 0) {
           tex = texture(textures[i], v_texcoord);
           mixTexture += tex;
      }
    }
    mixTexture = normalize(mixTexture);

    fragColor = clamp(mixTexture * 0.8 + mixTexture * 0.2 * NL, 0.0, 1.0);
}
