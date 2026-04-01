# -*- coding: utf-8 -*-

from abaqus import *
from abaqusConstants import *
import part
import mesh
import regionToolset
import numpy as np


def gerar_chassi_completo(L, h, b, r, R, E, t, peso, proporcao_banco):

    # ==========================================================
    # calculo das medidas
    # ==========================================================
    b_eff = b - (2 * r)

    h_eff = h
    if h_eff > 1e-6:
        h_eff -= (2*r)

    L_lat = L - (1.5 * R) - (2 * r)
    y_X = 0.4318181818181818 * L_lat
    reta_bico = 0.580574 * E
    bico_eff = (b_eff / 3.0) + r

    x_ext = b_eff / 2.0
    x_bico = bico_eff / 2.0

    comp_y_diagonal = L - (L_lat + reta_bico) - (2*r)
    # ==========================================================
    # definição dos pontos
    # ==========================================================
    pt_tras_esq = (-x_ext, 0.0, 0.0)
    pt_tras_dir = (x_ext, 0.0, 0.0)
    pt_lat_esq = (-x_ext, L_lat, 0.0)
    pt_lat_dir = (x_ext, L_lat, 0.0)
    pt_x_esq_inicio = (-x_ext, y_X, 0.0)
    pt_x_dir_inicio = (x_ext, y_X, 0.0)
    pt_alto_esq = (-x_ext, L_lat, h_eff)
    pt_alto_dir = (x_ext, L_lat, h_eff)
    y_reta = L_lat + reta_bico
    pt_reta_esq = (-x_ext, y_reta, h_eff)
    pt_reta_dir = (x_ext, y_reta, h_eff)
    y_frente = y_reta + comp_y_diagonal
    pt_bico_esq = (-x_bico, y_frente, h_eff)
    pt_bico_dir = (x_bico, y_frente, h_eff)

    # ==========================================================
    # 3. construção do wire
    # ==========================================================
    meu_modelo = mdb.models['Model-1']

    if 'chassi_parametrizado' in meu_modelo.parts:
        del meu_modelo.parts['chassi_parametrizado']

    minha_part = meu_modelo.Part(
        name='chassi_parametrizado', dimensionality=THREE_D, type=DEFORMABLE_BODY)

    # lista de barras
    linhas_wire = [
        (pt_tras_esq, pt_tras_dir),
        (pt_tras_esq, pt_lat_esq),
        (pt_tras_dir, pt_lat_dir),
        (pt_x_esq_inicio, pt_lat_dir),
        (pt_x_dir_inicio, pt_lat_esq),
        (pt_reta_esq, pt_bico_esq),
        (pt_reta_dir, pt_bico_dir),
        (pt_bico_esq, pt_bico_dir)
    ]

    # se a altura for maior q zero, desenha o degrau
    if h_eff > 1e-6:
        linhas_wire.extend([
            (pt_lat_esq, pt_alto_esq),
            (pt_lat_dir, pt_alto_dir),
            (pt_alto_esq, pt_reta_esq),
            (pt_alto_dir, pt_reta_dir)
        ])
    else:
        linhas_wire.extend([
            (pt_lat_esq, pt_reta_esq),
            (pt_lat_dir, pt_reta_dir)
        ])

    # Passa a lista dinâmica para o Abaqus desenhar
    minha_part.WirePolyLine(points=tuple(linhas_wire),
                            mergeType=IMPRINT, meshable=ON)

    # ==========================================================
    # cria um material pra conseguir criar a seção, e depois a seção
    # ==========================================================
    if 'PVC' not in meu_modelo.materials:
        meu_material = meu_modelo.Material(name='PVC')
        meu_material.Elastic(table=((3e9, 0.38), ))
        meu_material.Plastic(table=((52000000.0, 0.0), ))
        meu_material.Density(table=((1400.0, ), ))

    nome_perfil = 'Perfil_Tubo_R' + str(r)
    if nome_perfil not in meu_modelo.profiles:
        meu_modelo.PipeProfile(name=nome_perfil, r=r, t=t)

    if 'secao_chassi' not in meu_modelo.sections:
        meu_modelo.BeamSection(name='secao_chassi', integration=DURING_ANALYSIS,
                               poissonRatio=0.3, profile=nome_perfil,
                               material='PVC', temperatureVar=LINEAR)

    # ==========================================
    # orientação dos fios
    # ==========================================
    todas_arestas = minha_part.edges
    regiao_chassi = minha_part.Set(edges=todas_arestas, name='Todo_o_Chassi')
    minha_part.SectionAssignment(
        region=regiao_chassi, sectionName='secao_chassi')

    arestas_verticais, arestas_horizontais = [], []
    for aresta in todas_arestas:
        v1 = minha_part.vertices[aresta.getVertices()[0]].pointOn[0]
        v2 = minha_part.vertices[aresta.getVertices()[1]].pointOn[0]
        dz = abs(v1[2] - v2[2])
        dx = abs(v1[0] - v2[0])
        dy = abs(v1[1] - v2[1])

        if dz > (dx + dy):
            arestas_verticais.append(aresta)
        else:
            arestas_horizontais.append(aresta)

    if arestas_horizontais:
        reg_horiz = minha_part.Set(edges=part.EdgeArray(
            arestas_horizontais), name='Arestas_Horizontais')
        minha_part.assignBeamSectionOrientation(
            region=reg_horiz, method=N1_COSINES, n1=(0.0, 0.0, -1.0))

    if arestas_verticais:
        reg_vert = minha_part.Set(edges=part.EdgeArray(
            arestas_verticais), name='Arestas_Verticais')
        minha_part.assignBeamSectionOrientation(
            region=reg_vert, method=N1_COSINES, n1=(1.0, 0.0, 0.0))

    # ==========================================
    # 6. MALHA (B31)
    # ==========================================
    tipo_elemento = mesh.ElemType(elemCode=B31, elemLibrary=STANDARD)
    minha_part.setElementType(regions=regiao_chassi,
                              elemTypes=(tipo_elemento,))
    minha_part.seedPart(size=0.02, deviationFactor=0.1, minSizeFactor=0.1)
    minha_part.generateMesh()

    # ==========================================
    # 7. assembly, ponto de referência, acoplamento, carga e condições de contorno
    # ==========================================
    meu_assembly = meu_modelo.rootAssembly
    meu_assembly.DatumCsysByDefault(CARTESIAN)

    if 'Chassi' in meu_assembly.instances:
        del meu_assembly.instances['Chassi']

    instancia = meu_assembly.Instance(
        name='Chassi', part=minha_part, dependent=ON)

    # ponto da carga

    y_rp = L_lat * proporcao_banco
    coord_rp = (0.0, y_rp, h)

    rp_feature = meu_assembly.ReferencePoint(point=coord_rp)
    rp_id = rp_feature.id
    regiao_rp = meu_assembly.Set(referencePoints=(
        meu_assembly.referencePoints[rp_id],), name='Set_RP_Carga')

    # Busca a barra lateral em dois pontos (antes e depois do x cruzar a barra)
    # Ponto 1: Metade do caminho antes do X
    y_busca_1 = y_X / 2.0
    # Ponto 2: Metade do caminho entre o X e a dobra da haste vertical
    y_busca_2 = y_X + ((L_lat - y_X) / 2.0)

    # Encontra os pedaços da esquerda
    aresta_esq_1 = instancia.edges.findAt(((-x_ext, y_busca_1, 0.0),))
    aresta_esq_2 = instancia.edges.findAt(((-x_ext, y_busca_2, 0.0),))

    # Encontra os pedaços da direita
    aresta_dir_1 = instancia.edges.findAt(((x_ext, y_busca_1, 0.0),))
    aresta_dir_2 = instancia.edges.findAt(((x_ext, y_busca_2, 0.0),))

    # Soma todos os pedaços encontrados para formar as barras inteiras
    regiao_laterais = meu_assembly.Set(
        edges=aresta_esq_1 + aresta_esq_2 + aresta_dir_1 + aresta_dir_2,
        name='Set_Barras_Laterais')

    meu_modelo.Coupling(name='Acoplamento_Peso', controlPoint=regiao_rp, surface=regiao_laterais,
                        influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, weightingMethod=UNIFORM)

    if 'Step-Carga' not in meu_modelo.steps:
        meu_modelo.StaticStep(name='Step-Carga', previous='Initial', nlgeom=ON)

    peso_crianca = peso * -9.81
    meu_modelo.ConcentratedForce(
        name='Peso_Crianca', createStepName='Step-Carga', region=regiao_rp, cf3=peso_crianca)

    vertices_engaste = instancia.vertices.findAt(
        (pt_tras_esq,),
        (pt_tras_dir,),
        (pt_bico_esq,),
        (pt_bico_dir,)
    )

    regiao_apoios = meu_assembly.Set(vertices=vertices_engaste, name='apoios')
    meu_modelo.DisplacementBC(name='engaste', createStepName='Initial', region=regiao_apoios,
                              u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)

    print("Modelo 100% configurado! Geometria, RP, Carga e Malha gerados dinamicamente.")


# ==============================================================================================================================
# medidas físicas do chassi
# ==============================================================================================================================
L = 53.7         # distância entre eixos
h = 15   # altura da haste da roda
b = 38    # largura
R = 11    # raio da roda
E = 15.5        # distância entre o eixo frontal e a extremidade da roda mais próxima ao centro

r = 1.0     # raio da seção transversal
t = 0.15   # espessura da seção transversal (para o perfil tubular)

peso = 50  # peso da criança (em kg)

# distância (em porcentagem) do eixo traseiro até o ponto onde o peso é aplicado na barra lateral
proporcao_banco = 0.5

variaveis = np.array([L, h, b, r, R, E, t])
variaveis *= 0.01  # Converte para metros

gerar_chassi_completo(*variaveis, peso=peso, proporcao_banco=proporcao_banco)
