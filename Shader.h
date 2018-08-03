#pragma once

#include <cstdint>

class Shader
{
public:
    enum Signature {
        ATTRIBUTE_LESS = -1,
        STANDARD = 0,
        COMPATIBLE = 1,
        MAX
    };

    Shader() {}
    ~Shader() {}

    bool loadVertexShader(const char* filename);
    bool loadFragmentShader(const char* filename);
    bool compileVertexShader(char* buffer);
    bool compileFragmentShader(char* buffer);
    void compileShader(uint32_t type, const char* code);

    bool createProgram(int signature = STANDARD);
    void destroyProgram();

    inline uint32_t getProgram() { return m_ProgramObject; }

    static void SetVersionString(int32_t version);
    
public:
    static char g_VersionString[16];
    static uint32_t g_VersionStringSize;
    
    static char* LoadFileToMemory(const char* filename);
   
private:
    uint32_t m_ProgramObject;
    uint32_t m_VertexShader;
    uint32_t m_FragmentShader;
};

