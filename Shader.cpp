
#include "Shader.h"

// TODO: replace by (f)iostream
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

// TODO: common file
#if defined(_WIN32)
#include <GL/glew.h>
#include <GL/GL.h>
#elif defined(__APPLE__)
#include <OpenGL/OpenGL.h>
#include <OpenGL/gl3.h>
#include <OpenGL/gl3ext.h>
#endif

// ---

struct Include
{
    char* buffer;
    int len;
    char filename[100];
    
    static std::vector<Include> includeFiles;
    
    static int Exists(char* fileName) {
        int index = 0;
        for (auto& i : includeFiles) {
            if (strcmp(i.filename, fileName) == 0) {
                return index;
            }
            index++;
        }
        return -1;
    }
};
std::vector<Include> Include::includeFiles;

char Shader::g_VersionString[16];
uint32_t Shader::g_VersionStringSize = 0;

void Shader::SetVersionString(int32_t version)
{
#if defined(TARGET_IOS) || defined(TARGET_ANDROID)
    snprintf(Shader::g_VersionString, 16, "#version %d es\n", version);
#else
    snprintf(Shader::g_VersionString, 16, "#version %d\n", version);
#endif
    std::cout << Shader::g_VersionString << std::endl;
    Shader::g_VersionStringSize = (uint32_t)strnlen(Shader::g_VersionString, 16);
}

// ---

void Shader::destroyProgram()
{
    glDetachShader(m_ProgramObject, m_VertexShader);
    glDetachShader(m_ProgramObject, m_FragmentShader);
    glDeleteProgram(m_ProgramObject);
    glDeleteShader(m_VertexShader);
    glDeleteShader(m_FragmentShader);
}

bool CheckProgram(GLuint program)
{
    int logLength = 0;
    int didLink = 0;
    glGetProgramiv(program, GL_LINK_STATUS, &didLink);
    if (!didLink)
    {
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLength);
        if (logLength) {
            char logBuffer[4096];
            int len = 0;
            glGetProgramInfoLog(program, 4096, &len, logBuffer);
            std::cout << logBuffer << std::endl;
        }
    }
    return didLink == 1;
}

bool Shader::createProgram(int signature)
{
    m_ProgramObject = glCreateProgram();
    glAttachShader(m_ProgramObject, m_VertexShader);
    
    // TEMP: prefer layout(location) on ES3.0 and GL3.2+ (if ARB_explicit_attrib_location is available)
    if (signature == COMPATIBLE) {
        glBindAttribLocation(m_ProgramObject, 0, "a_Position");
        //if (glGetAttribLocation(m_ProgramObject, "a_Normal") >= 0)
            glBindAttribLocation(m_ProgramObject, 1, "a_Normal");
        //if (glGetAttribLocation(m_ProgramObject, "a_Color") >= 0)
            glBindAttribLocation(m_ProgramObject, 2, "a_Color");
        //if (glGetAttribLocation(m_ProgramObject, "a_TexCoords") >= 0)
            glBindAttribLocation(m_ProgramObject, 3, "a_TexCoords");
        //if (glGetAttribLocation(m_ProgramObject, "a_TexCoords1") >= 0)
        //    glBindAttribLocation(m_ProgramObject, 4, "a_TexCoords1");
    }
    glAttachShader(m_ProgramObject, m_FragmentShader);
    glLinkProgram(m_ProgramObject);
    CheckProgram(m_ProgramObject);
    return true;
}

// ---

// debug crap
static const char* processedShaderFileName;

#ifdef TARGET_ANDROID
char* LoadFileToMemory(const char* filename, int& len)
{
    extern int androidLoadFileToMemory(const char* filePath, char** buffer, const char* prefix = NULL, int prefixLen = 0);
    char* buffer = NULL;
    len = androidLoadFileToMemory(filename, &buffer);//, Shader::g_VersionString, Shader::g_VersionStringSize);
    return buffer;
}
#else
char* LoadFileToMemory(const char* filename, int& len)
{
    FILE* file = fopen(filename, "rb");
    fseek(file, 0, SEEK_END);
    len = (int)ftell(file);
    auto bufferLen = len;
    rewind(file);
    auto buffer = new char[bufferLen+1];
    fread(buffer, len, 1, file);
    buffer[bufferLen] = '\0'; // = 0
    fclose(file);
    
    return buffer;
}
#endif

#ifdef TARGET_ANDROID
char* Shader::LoadFileToMemory(const char* filename)
{
    extern int androidLoadFileToMemory(const char* filePath, char** buffer, const char* prefix = NULL, int prefixLen = 0);
    char* buffer = NULL;
    androidLoadFileToMemory(filename, &buffer, g_VersionString, g_VersionStringSize);
    return buffer;
}
#else
char* Shader::LoadFileToMemory(const char* filename)
{
    FILE* file = fopen(filename, "rb");
    fseek(file, 0, SEEK_END);
    auto len = ftell(file);
    auto bufferLen = g_VersionStringSize+len;
    rewind(file);
    auto buffer = new char[bufferLen+1];
    strncpy(buffer, g_VersionString, g_VersionStringSize);
    fread(buffer+g_VersionStringSize, len, 1, file);
    buffer[bufferLen] = '\0'; // = 0
    fclose(file);
    
    return buffer;
}
#endif

bool CheckShader(GLuint shader)
{
    int logLength = 0;
    int didCompile = 0;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLength);
    glGetShaderiv(shader, GL_COMPILE_STATUS, &didCompile);
    if (!didCompile || logLength)
    {
        if (processedShaderFileName)
            std::cout << (processedShaderFileName) << std::endl;
        char logBuffer[4096];
        int len = 0;
        glGetShaderInfoLog(shader, 4096, &len, logBuffer);
        std::cout << (logBuffer) << std::endl;
    }
    return didCompile == 1;
}

void Shader::compileShader(uint32_t type, const char* code)
{
    uint32_t shaderObject = glCreateShader(type);
    if (shaderObject != 0) {
        glShaderSource(shaderObject, 1, &code, nullptr);
        glCompileShader(shaderObject);
        CheckShader(shaderObject);
        switch (type)
        {
            case GL_VERTEX_SHADER: m_VertexShader = shaderObject; break;
            case GL_FRAGMENT_SHADER: m_FragmentShader = shaderObject; break;
            default: break;
        }
    }
}

char* skipWhiteSpace(char* cursor)
{
    while (*cursor == ' ' || *cursor == '\t')
        cursor++;
    return (cursor);
}

char* skipLine(char* cursor)
{
    while (*cursor != '\n' && *cursor != 0)
        cursor++;
    return (++cursor);
}

char* findDoubleQuote(char* cursor)
{
    while (*cursor != '\"')
        cursor++;
    return (++cursor);
}

char* skipDirective(char* cursor)
{
    char* oldCursor = cursor;
    //while (*cursor != 0)
    {
        char *currentLine = cursor;
        cursor = skipWhiteSpace(cursor);
        if (*cursor == '#') {
            return skipLine(cursor);
        }
    }
    return oldCursor;
}

char* findInclude(char* cursor)
{
    while (*cursor != 0)
    {
        char *currentLine = cursor;
        cursor = skipWhiteSpace(cursor);
        if (*cursor == '#') {
            if (strncmp(cursor, "#include", 8) == 0) {
                return currentLine;
            }
            cursor = skipDirective(cursor);
        }
        else if ((*cursor == '/') && (*(cursor+1) == '/')) {
            //skip comments
            cursor = skipLine(cursor);
        }
        else {
            cursor = skipLine(cursor);
        }
    }
    return nullptr;
}

// ---

bool resolveIncludes(char** buffer)
{
    struct Occurence {
        char* cursor;
        int id;
    };
    Occurence includeOccurence[32];
    int occurences = 0;
    
    bool hasInclude = false;
    
    // Look for include references & preload files
    char* cursor = *buffer;
    do {
        cursor = findInclude(cursor);
        if (cursor) {
            char* includeLine = cursor;
            char* filename = findDoubleQuote(cursor);
            cursor = findDoubleQuote(filename);
            char* end = cursor-1;
            //*(end) = 0; // marker
            int index = Include::Exists(filename);
            if (index == -1) {
                Include include;
                strncpy(include.filename, filename, end-filename);
                include.filename[end-filename] = 0;
                char filepath[1024];
                strcpy(filepath, ""/*Engine::PathToResources*/);
                strcat(filepath, include.filename);
                include.buffer = LoadFileToMemory(filepath, include.len);
                index = (int)Include::includeFiles.size();
                Include::includeFiles.push_back(include);
            }
            
            includeOccurence[occurences].cursor = includeLine;
            includeOccurence[occurences].id = index;
            occurences++;
            
            //cursor = skipLine(cursor);
            hasInclude = true;
        }
    }
    while (cursor);
    
    if (hasInclude)
    {
        char* oldBuffer = *buffer;
        
        cursor = oldBuffer;
        
        int includeSizes = 0;
        for (int i = 0; i < occurences; i++) {
            includeSizes += Include::includeFiles[includeOccurence[i].id].len;
        }
        // todo, calc strnlen during first pass
        int oldBufferLen = (int)strlen(*buffer);
        char* output = new char[oldBufferLen + includeSizes];
        *buffer = output;
        
        char* prevPos = cursor;
        for (int i = 0; i < occurences; i++) {
            // copy unchanged text before include directive
            char* incPos = includeOccurence[i].cursor;
            int pretextLen = (int)(incPos - prevPos);
            if (pretextLen > 0) {
                memcpy(output, prevPos, pretextLen);
            }
            prevPos = skipDirective(incPos);
            cursor = prevPos;
            output += pretextLen;
            // replace include directive
            
            Include& include = Include::includeFiles[includeOccurence[i].id];
            memcpy(output, include.buffer, include.len);
            output += include.len;
        }
        // last part
        cursor = skipDirective(cursor);
        uint32_t remaining = (uint32_t)((oldBuffer+oldBufferLen) - cursor + 1);
        memcpy(output, cursor, remaining);
        
        delete[] oldBuffer;
        
    }
    
    return hasInclude;
}

bool Shader::loadVertexShader(const char* filename)
{
    processedShaderFileName = filename;
    char filepath[1024];
    strcpy(filepath, ""/*Engine::PathToResources*/);
    strcat(filepath, filename);
    auto buffer = LoadFileToMemory(filepath);
    bool has = resolveIncludes(&buffer);
    compileShader(GL_VERTEX_SHADER, buffer);
    delete[] buffer;
    
    return (m_VertexShader != 0);
}

bool Shader::loadFragmentShader(const char* filename)
{
    processedShaderFileName = filename;
    char filepath[1024];
    strcpy(filepath, ""/*Engine::PathToResources*/);
    strcat(filepath, filename);
    auto buffer = LoadFileToMemory(filepath);
    bool has = resolveIncludes(&buffer);
    compileShader(GL_FRAGMENT_SHADER, buffer);
    delete[] buffer;

    return (m_FragmentShader != 0);
}

bool Shader::compileVertexShader(char* buffer)
{
    bool has = resolveIncludes(&buffer);
    compileShader(GL_VERTEX_SHADER, buffer);

    return (m_VertexShader != 0);
}

bool Shader::compileFragmentShader(char* buffer)
{
    bool has = resolveIncludes(&buffer);
    compileShader(GL_FRAGMENT_SHADER, buffer);

    return (m_FragmentShader != 0);
}

