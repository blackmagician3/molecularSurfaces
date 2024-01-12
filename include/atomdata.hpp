#ifndef ATOMDATA_HPP
#define ATOMDATA_HPP

#include <string>

#ifndef ELEMENTS
#define ELEMENTS 112
#endif

typedef unsigned int uint;

struct entry
{
    // element symbol
    std::string symbol;
    // radii
    float rSingleBonds1;
    float rSingleBonds2;
    float rVanDerWaals;
    // colors depending on color scheme
    uint cCorey;
    uint cKoltun;
    uint cJmol;
    uint cRasmolOld;
    uint cRasmolNew;
    uint cPubChem;
};

const entry pse[ELEMENTS] = {
    {"H", .31f, .32f, 1.2f, (uint)0xffffff, (uint)0xffffff, (uint)0xffffff, (uint)0xffffff, (uint)0xffffff, (uint)0x638c8c},
    {"He", .28f, .46f, 1.4f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xd9ffff, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xd593a1},
    {"Li", 1.28f, 1.33f, 1.82f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xcc80ff, (uint)0xb22222, (uint)0xb22121, (uint)0xd56632},
    {"Be", .96f, 1.02f, 1.53f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xc2ff00, (uint)0xff1493, (uint)0xff1493, (uint)0xd5bad5},
    {"B", .84f, .85f, 1.92f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffb5b5, (uint)0x00ff00, (uint)0x00ff00, (uint)0x2ad52a},
    {"C", .76f, .75f, 1.7f, (uint)0x202020, (uint)0x202020, (uint)0x909090, (uint)0xc8c8c8, (uint)0xd3d3d3, (uint)0x274a4a},
    {"N", .71f, .71f, 1.55f, (uint)0x2060ff, (uint)0x2060ff, (uint)0x3050f8, (uint)0x8f8fff, (uint)0x87cee6, (uint)0x0000ff},
    {"O", .66f, .63f, 1.52f, (uint)0xee2010, (uint)0xee2010, (uint)0xff0d0d, (uint)0xf00000, (uint)0xff0000, (uint)0xff0000},
    {"F", .57f, .64f, 1.35f, (uint)0xffc0cb, (uint)0x00ff00, (uint)0x90e050, (uint)0xdaa520, (uint)0xdaa520, (uint)0xd52092},
    {"Ne", .58f, .67f, 1.54f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xb3e3f5, (uint)0xff1493, (uint)0xff1493, (uint)0xff00ff},
    {"Na", 1.66f, 1.55f, 2.27f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xab5cf2, (uint)0x0000ff, (uint)0x0000ff, (uint)0x0e73d5},
    {"Mg", 1.41f, 1.39f, 1.73f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x8aff00, (uint)0x228b22, (uint)0x228b22, (uint)0x198c19},
    {"Al", 1.21f, 1.26f, 1.84f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xbfa6a6, (uint)0x808090, (uint)0x696969, (uint)0x838c8c},
    {"Si", 1.11f, 1.16f, 2.1f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xf0c8a0, (uint)0xdaa520, (uint)0xdaa520, (uint)0xd59e13},
    {"P", 1.07f, 1.11f, 1.8f, (uint)0xffc0cb, (uint)0x8020ff, (uint)0xff8000, (uint)0xffa500, (uint)0xffaa00, (uint)0xd58600},
    {"S", 1.05f, 1.03f, 1.8f, (uint)0xffc0cb, (uint)0xffff00, (uint)0xffff30, (uint)0xffc832, (uint)0xffff00, (uint)0xd5d500},
    {"Cl", 1.02f, .99f, 1.75f, (uint)0xffc0cb, (uint)0x00bb00, (uint)0x1ff01f, (uint)0x00ff00, (uint)0x00ff00, (uint)0x2ad52a},
    {"Ar", 1.06f, .96f, 1.88f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x80d1e3, (uint)0xff1493, (uint)0xff1493, (uint)0xff00ff},
    {"K", 2.03f, 1.96f, 2.75f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x8f40d4, (uint)0xff1493, (uint)0xff1493, (uint)0xd50575},
    {"Ca", 1.76f, 1.71f, 2.31f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x3dff00, (uint)0x808090, (uint)0x696969, (uint)0x838c8c},
    {"Sc", 1.7f, 1.48f, 2.11f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xe6e6e6, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Ti", 1.6f, 1.36f, 1.87f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xbfc2c7, (uint)0x808090, (uint)0x696969, (uint)0x838c8c},
    {"v", 1.53f, 1.34f, 1.79f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xa6a6ab, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Cr", 1.39f, 1.22f, 1.89f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x8a99c7, (uint)0x808090, (uint)0x696969, (uint)0x838c8c},
    {"Mn", 1.5f, 1.19f, 1.97f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x9c7ac7, (uint)0x808090, (uint)0x696969, (uint)0x838c8c},
    {"Fe", 1.42f, 1.16f, 1.94f, (uint)0xffc0cb, (uint)0xd0d0d0, (uint)0xe06633, (uint)0xffa500, (uint)0xffaa00, (uint)0xffa900},
    {"Co", 1.38f, 1.11f, 1.92f, (uint)0xffc0cb, (uint)0xd0d0d0, (uint)0xf090a0, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Ni", 1.24f, 1.1f, 1.63f, (uint)0xffc0cb, (uint)0xd0d0d0, (uint)0x50d050, (uint)0xa52a2a, (uint)0x802828, (uint)0x838c8c},
    {"Cu", 1.32f, 1.12f, 1.4f, (uint)0xffc0cb, (uint)0xd0d0d0, (uint)0xc88033, (uint)0xa52a2a, (uint)0x802828, (uint)0x838c8c},
    {"Zn", 1.22f, 1.18f, 1.39f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x7d80b0, (uint)0xa52a2a, (uint)0x802828, (uint)0x838c8c},
    {"Ga", 1.22f, 1.24f, 1.87f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xc28f8f, (uint)0xff1493, (uint)0xff1493, (uint)0xd5cd72},
    {"Ge", 1.2f, 1.21f, 2.11f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x668f8f, (uint)0xff1493, (uint)0xff1493, (uint)0xd5cd72},
    {"As", 1.19f, 1.21f, 1.85f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xbd80e3, (uint)0xff1493, (uint)0xff1493, (uint)0xd56632},
    {"Se", 1.2f, 1.16f, 1.9f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffa100, (uint)0xff1493, (uint)0xff1493, (uint)0xd5cd72},
    {"Br", 1.2f, 1.14f, 1.83f, (uint)0xffc0cb, (uint)0x008800, (uint)0xa62929, (uint)0xa52a2a, (uint)0x802828, (uint)0xd58639},
    {"Kr", 1.16f, 1.17f, 2.02f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x5cb8d1, (uint)0xff1493, (uint)0xff1493, (uint)0xff00ff},
    {"Rb", 2.2f, 2.1f, 3.03f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x702eb0, (uint)0xff1493, (uint)0xff1493, (uint)0xff00ff},
    {"Sr", 1.95f, 1.85f, 2.49f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x00ff00, (uint)0xff1493, (uint)0xff1493, (uint)0xff0000},
    {"Y", 1.9f, 1.63f, 2.19f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x94ffff, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Zr", 1.75f, 1.54f, 1.86f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x94e0e0, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Nb", 1.64f, 1.47f, 2.07f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x73c2c9, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Mo", 1.54f, 1.38f, 2.09f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x54b5b5, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Tc", 1.47f, 1.28f, 2.09f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x3b9e9e, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Ru", 1.46f, 1.25f, 2.07f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x248f8f, (uint)0xff1493, (uint)0xff1493, (uint)0xd5b100},
    {"Rh", 1.42f, 1.25f, 1.95f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x0a7d8c, (uint)0xff1493, (uint)0xff1493, (uint)0xd5b100},
    {"Pd", 1.39f, 1.2f, 2.02f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x006985, (uint)0xff1493, (uint)0xff1493, (uint)0xd5b100},
    {"Ag", 1.45f, 1.28f, 1.72f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xc0c0c0, (uint)0x808090, (uint)0x696969, (uint)0x838c8c},
    {"Cd", 1.44f, 1.36f, 1.58f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffd98f, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"In", 1.42f, 1.42f, 1.93f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xa67573, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Sn", 1.39f, 1.4f, 2.17f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x668080, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Sb", 1.39f, 1.4f, 2.06f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x9e63b5, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Te", 1.38f, 1.36f, 2.06f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xd47a00, (uint)0xff1493, (uint)0xff1493, (uint)0xd5cd72},
    {"I", 1.39f, 1.33f, 1.98f, (uint)0xffc0cb, (uint)0x005500, (uint)0x940094, (uint)0xa020f0, (uint)0xa020f0, (uint)0xff00ff},
    {"Xe", 1.4f, 1.31f, 2.16f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x429eb0, (uint)0xff1493, (uint)0xff1493, (uint)0xff00ff},
    {"Cs", 2.44f, 2.32f, 3.43f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x57178f, (uint)0xff1493, (uint)0xff1493, (uint)0xd598d5},
    {"Ba", 2.15f, 1.96f, 2.68f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x00c900, (uint)0xffa500, (uint)0xffaa00, (uint)0x2ad52a},
    {"La", 2.07f, 1.8f, 2.4f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x70d4ff, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Ce", 2.04f, 1.63f, 2.35f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffffc7, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Pr", 2.03f, 1.76f, 2.39f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xd9ffc7, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Nd", 2.01f, 1.74f, 2.29f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xc7ffc7, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Pm", 1.99f, 1.73f, 2.36f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xa3ffc7, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Sm", 1.98f, 1.72f, 2.29f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x8fffc7, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Eu", 1.98f, 1.68f, 2.33f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x61ffc7, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Gd", 1.96f, 1.69f, 2.37f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x45ffc7, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Tb", 1.94f, 1.68f, 2.21f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x30ffc7, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Dy", 1.92f, 1.67f, 2.29f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x1fffc7, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Ho", 1.92f, 1.66f, 2.16f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x00ff9c, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Er", 1.89f, 1.65f, 2.35f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x00e675, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Tm", 1.9f, 1.64f, 2.27f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x00d452, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Yb", 1.87f, 1.7f, 2.42f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x00bf38, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Lu", 1.87f, 1.62f, 2.21f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x00ab24, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Hf", 1.75f, 1.52f, 2.12f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x4dc2ff, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Ta", 1.7f, 1.46f, 2.17f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x4da6ff, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"W", 1.62f, 1.37f, 2.1f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x2194d6, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Re", 1.51f, 1.31f, 2.17f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x267dab, (uint)0xff1493, (uint)0xff1493, (uint)0xd5b100},
    {"Os", 1.44f, 1.29f, 2.16f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x266696, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Ir", 1.41f, 1.22f, 2.02f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x175487, (uint)0xff1493, (uint)0xff1493, (uint)0xd5b100},
    {"Pt", 1.36f, 1.23f, 2.09f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xd0d0e0, (uint)0xff1493, (uint)0xff1493, (uint)0xd5b100},
    {"Au", 1.36f, 1.24f, 1.66f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffd123, (uint)0xdaa520, (uint)0xdaa520, (uint)0xd5b100},
    {"Hg", 1.32f, 1.33f, 2.09f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xb8b8d0, (uint)0xff1493, (uint)0xff1493, (uint)0x6e8092},
    {"Tl", 1.45f, 1.44f, 1.96f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xa6544d, (uint)0xff1493, (uint)0xff1493, (uint)0xd593a1},
    {"Pb", 1.46f, 1.44f, 2.02f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x575961, (uint)0xff1493, (uint)0xff1493, (uint)0x838c8c},
    {"Bi", 1.48f, 1.51f, 2.07f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x9e4fb5, (uint)0xff1493, (uint)0xff1493, (uint)0xd593a1},
    {"Po", 1.4f, 1.45f, 1.97f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xab5c00, (uint)0xff1493, (uint)0xff1493, (uint)0xd593a1},
    {"At", 1.5f, 1.47f, 2.02f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x754f45, (uint)0xff1493, (uint)0xff1493, (uint)0xff00ff},
    {"Rn", 1.5f, 1.42f, 2.2f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x428296, (uint)0xff1493, (uint)0xff1493, (uint)0xd598d5},
    {"Fr", 2.6f, 2.23f, 3.48f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x420066, (uint)0xff1493, (uint)0xff1493, (uint)0xd598d5},
    {"Ra", 2.21f, 2.01f, 2.83f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x007d00, (uint)0xff1493, (uint)0xff1493, (uint)0x2ad52a},
    {"Ac", 2.15f, 1.86f, 2.6f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x70abfa, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Th", 2.06f, 1.75f, 2.37f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x00baff, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Pa", 2.0f, 1.69f, 2.43f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x00a1ff, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"U", 1.96f, 1.7f, 2.4f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x008fff, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Np", 1.9f, 1.71f, 2.21f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x0080ff, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Pu", 1.87f, 1.72f, 2.43f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x006bff, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Am", 1.8f, 1.66f, 2.44f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x545cf2, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Cm", 1.69f, 1.66f, 2.45f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x785ce3, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Bk", .0f, 1.68f, 2.44f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0x8a4fe3, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Cf", .0f, 1.68f, 2.45f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xa136d4, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Es", .0f, 1.65f, 2.45f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xb31fd4, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Fm", .0f, 1.67f, .0f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xb31fba, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Md", .0f, 1.73f, .0f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xb30da6, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"No", .0f, 1.76f, .0f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xbd0d87, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Lr", .0f, 1.61f, .0f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xc70066, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Rf", .0f, 1.57f, .0f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xcc0059, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Db", .0f, 1.49f, .0f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xd1004f, (uint)0xff1493, (uint)0xff1493, (uint)0x00ccd5},
    {"Sg", .0f, 1.43f, .0f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffb6c1, (uint)0xffb6c1, (uint)0xffb6c1},
    {"Bh", .0f, 1.41f, .0f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffb6c1, (uint)0xffb6c1, (uint)0xffb6c1},
    {"Hs", .0f, 1.34f, .0f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffb6c1, (uint)0xffb6c1, (uint)0xffb6c1},
    {"Mt", .0f, 1.29f, .0f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffb6c1, (uint)0xffb6c1, (uint)0xffb6c1},
    {"Ds", .0f, 1.28f, .0f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffb6c1, (uint)0xffb6c1, (uint)0xffb6c1},
    {"Rg", .0f, 1.21f, .0f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffb6c1, (uint)0xffb6c1, (uint)0xffb6c1},
    {"Cn", .0f, 1.22f, .0f, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffc0cb, (uint)0xffb6c1, (uint)0xffb6c1, (uint)0xffb6c1}};

void findEntry(std::string symbol, uint colorScheme, float *radius_out, uint *color_out)
{

    uint n = 0;
    for (uint i = 0; i < ELEMENTS; i++)
    {
        // compare symbols
        if (symbol.compare(pse[i].symbol) == 0)
        {

            n = i;
            break;
        }
    }

    *radius_out = pse[n].rVanDerWaals;

    switch (colorScheme)
    {
    case 1:
        *color_out = pse[n].cCorey;
        break;
    case 2:
        *color_out = pse[n].cKoltun;
        break;
    case 3:
        *color_out = pse[n].cJmol;
        break;
    case 4:
        *color_out = pse[n].cRasmolOld;
        break;
    case 5:
        *color_out = pse[n].cRasmolNew;
        break;
    case 6:
        *color_out = pse[n].cPubChem;
        break;
    default:
        *color_out = 0x000000;
        break;
    }
}
#endif