# -*- coding: utf-8 -*-

from abaqus import *
from abaqusConstants import *
import part
import mesh
import regionToolset
import numpy as np


class GerarModeloChassi():
    def __init__(self, L, h, b, D, t, R, peso, proporcao_banco, nome, angulo_curva_bico=45, comp_joelho_cm=3.5):
        k = 0.01
        self.L = L*k
        self.h = h*k
        self.b = b*k
        self.D = D*k
        self.t = t*k

        # R representa o raio da roda do carrinho
        self.R = R*k

        # A medida E substitui a antiga constante empírica, utilizando L, R e a folga fixa de 1 cm
        self.folga = 1.0 * k
        self.E = self.L - self.R - self.folga

        self.peso = peso
        self.proporcao_banco = proporcao_banco
        self.angulo_bico = angulo_curva_bico

        # Fator de compensação das conexões (joelhos) de PVC
        self.comp_joelho = comp_joelho_cm * k

        if nome not in mdb.models:
            self.meu_modelo = mdb.Model(name=nome)
        else:
            self.meu_modelo = mdb.models[nome]
       

    def criar_material_base(self):
        if 'PVC' not in self.meu_modelo.materials:
            meu_material = self.meu_modelo.Material(name='PVC')
            meu_material.Elastic(table=((3e9, 0.38), ))
            meu_material.Plastic(table=((52000000.0, 0.0), ))
            meu_material.Density(table=((1400.0, ), ))

        nome_perfil = 'Perfil_Tubo_R' + str(self.D)
        if nome_perfil not in self.meu_modelo.profiles:
            self.meu_modelo.PipeProfile(name=nome_perfil, r=self.D, t=self.t)

        if 'secao_chassi' not in self.meu_modelo.sections:
            self.meu_modelo.BeamSection(name='secao_chassi', integration=DURING_ANALYSIS,
                                        poissonRatio=0.3, profile=nome_perfil,
                                        material='PVC', temperatureVar=LINEAR)

    def geometria(self):

        # Medidas do chassi com as compensações necessárias

        b_eff = self.b - (1 * self.D) - (2 * self.comp_joelho)     # largura
        h_eff = self.h                                             # altura da haste

        if h_eff > 1e-6:
            # faz a compensação de raio e conexões depois de verificar que h != 0
            h_eff -= (1 * self.D) + (2 * self.comp_joelho)

        # Comprimento efetivo descontando os joelhos da frente e de trás
        L_eff = self.L - (2 * self.comp_joelho)

        self.L_lat = L_eff - self.R - (1 * self.D)      # longarinas

        # 4.1818... foi constante encontrada nos cálculos (pode ser alterada para deslocar o X)
        self.y_X = 0.4318181818181818 * self.L_lat

        # reta lateral do bico utilizando o novo E calculado
        reta_bico = 0.580574 * self.R

        # largura da base, simplesmente pra facilitar os calculos
        self.x_ext = b_eff / 2.0
        comp_y_diagonal = L_eff - (self.L_lat + reta_bico) - (1*self.D)

        delta_x = comp_y_diagonal / np.tan(np.deg2rad(self.angulo_bico))

        # O x do bico será a borda externa menos a variação lateral
        x_bico = self.x_ext - delta_x

        # ==========================================================
        # definição dos pontos
        # ==========================================================
        y_reta = self.L_lat + reta_bico
        y_frente = y_reta + comp_y_diagonal

        # serão utilizados para engaste
        self.pt_tras_esq = (-self.x_ext, 0.0, 0.0)
        self.pt_tras_dir = (self.x_ext, 0.0, 0.0)
        self.pt_bico_esq = (-x_bico, y_frente, h_eff)
        self.pt_bico_dir = (x_bico, y_frente, h_eff)

        pt_lat_esq = (-self.x_ext, self.L_lat, 0.0)
        pt_lat_dir = (self.x_ext, self.L_lat, 0.0)
        pt_x_esq_inicio = (-self.x_ext, self.y_X, 0.0)
        pt_x_dir_inicio = (self.x_ext, self.y_X, 0.0)
        pt_alto_esq = (-self.x_ext, self.L_lat, h_eff)
        pt_alto_dir = (self.x_ext, self.L_lat, h_eff)

        pt_reta_esq = (-self.x_ext, y_reta, h_eff)
        pt_reta_dir = (self.x_ext, y_reta, h_eff)

        # ================================================================
        # Construção do wire

        if 'chassi_parametrizado' in self.meu_modelo.parts:
            del self.meu_modelo.parts['chassi_parametrizado']

        self.minha_part = self.meu_modelo.Part(
            name='chassi_parametrizado', dimensionality=THREE_D, type=DEFORMABLE_BODY)

        # lista de barras
        linhas_wire = [
            (self.pt_tras_esq, self.pt_tras_dir),
            (self.pt_tras_esq, pt_lat_esq),
            (self.pt_tras_dir, pt_lat_dir),
            (pt_x_esq_inicio, pt_lat_dir),
            (pt_x_dir_inicio, pt_lat_esq),
            (pt_reta_esq, self.pt_bico_esq),
            (pt_reta_dir, self.pt_bico_dir),
            (self.pt_bico_esq, self.pt_bico_dir)
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

        # desenha o wire
        self.minha_part.WirePolyLine(points=tuple(
            linhas_wire),  mergeType=IMPRINT, meshable=ON)

        # ==========================================
        # orientação dos fios
        # ==========================================

        self.criar_material_base()

        todas_arestas = self.minha_part.edges

        self.regiao_chassi = self.minha_part.Set(
            edges=todas_arestas, name='Todo_o_Chassi')
        self.minha_part.SectionAssignment(
            region=self.regiao_chassi, sectionName='secao_chassi')

        arestas_verticais, arestas_horizontais = [], []
        for aresta in todas_arestas:
            v1 = self.minha_part.vertices[aresta.getVertices()[0]].pointOn[0]
            v2 = self.minha_part.vertices[aresta.getVertices()[1]].pointOn[0]
            dz = abs(v1[2] - v2[2])
            dx = abs(v1[0] - v2[0])
            dy = abs(v1[1] - v2[1])

            if dz > (dx + dy):
                arestas_verticais.append(aresta)
            else:
                arestas_horizontais.append(aresta)

        if arestas_horizontais:
            reg_horiz = self.minha_part.Set(edges=part.EdgeArray(
                arestas_horizontais), name='Arestas_Horizontais')
            self.minha_part.assignBeamSectionOrientation(
                region=reg_horiz, method=N1_COSINES, n1=(0.0, 0.0, -1.0))

        if arestas_verticais:
            reg_vert = self.minha_part.Set(edges=part.EdgeArray(
                arestas_verticais), name='Arestas_Verticais')
            self.minha_part.assignBeamSectionOrientation(
                region=reg_vert, method=N1_COSINES, n1=(1.0, 0.0, 0.0))

    def gerar_malha(self):
        tipo_elemento = mesh.ElemType(elemCode=B31, elemLibrary=STANDARD)
        self.minha_part.setElementType(
            regions=self.regiao_chassi, elemTypes=(tipo_elemento,))
        self.minha_part.seedPart(
            size=0.02, deviationFactor=0.1, minSizeFactor=0.1)
        self.minha_part.generateMesh()

    def assembly(self):
        self.meu_assembly = self.meu_modelo.rootAssembly
        self.meu_assembly.DatumCsysByDefault(CARTESIAN)

        if 'Chassi' in self.meu_assembly.instances:
            del self.meu_assembly.instances['Chassi']

        self.instancia = self.meu_assembly.Instance(
            name='Chassi', part=self.minha_part, dependent=ON)

    def carga_engaste(self):
        y_rp = self.L_lat * self.proporcao_banco
        coord_rp = (0.0, y_rp, self.h)

        rp_feature = self.meu_assembly.ReferencePoint(point=coord_rp)
        rp_id = rp_feature.id
        regiao_rp = self.meu_assembly.Set(referencePoints=(
            self.meu_assembly.referencePoints[rp_id],), name='Set_RP_Carga')

        # Busca a barra lateral em dois pontos (antes e depois do x cruzar a barra)
        # Ponto 1: Metade do caminho antes do X
        y_busca_1 = self.y_X / 2.0
        # Ponto 2: Metade do caminho entre o X e a dobra da haste vertical
        y_busca_2 = self.y_X + ((self.L_lat - self.y_X) / 2.0)

        # Encontra os pedaços da esquerda
        aresta_esq_1 = self.instancia.edges.findAt(
            ((-self.x_ext, y_busca_1, 0.0),))
        aresta_esq_2 = self.instancia.edges.findAt(
            ((-self.x_ext, y_busca_2, 0.0),))

        # Encontra os pedaços da direita
        aresta_dir_1 = self.instancia.edges.findAt(
            ((self.x_ext, y_busca_1, 0.0),))
        aresta_dir_2 = self.instancia.edges.findAt(
            ((self.x_ext, y_busca_2, 0.0),))

        # Soma todos os pedaços encontrados para formar as barras inteiras
        regiao_laterais = self.meu_assembly.Set(
            edges=aresta_esq_1 + aresta_esq_2 + aresta_dir_1 + aresta_dir_2,
            name='Set_Barras_Laterais')

        self.meu_modelo.Coupling(name='Acoplamento_Peso', controlPoint=regiao_rp, surface=regiao_laterais,
                                 influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, weightingMethod=UNIFORM)

        if 'Step-Carga' not in self.meu_modelo.steps:
            self.meu_modelo.StaticStep(
                name='Step-Carga', previous='Initial', nlgeom=ON)

        peso_crianca = self.peso * -9.81

        self.meu_modelo.ConcentratedForce(
            name='Peso_Crianca', createStepName='Step-Carga', region=regiao_rp, cf3=peso_crianca)

        vertices_engaste = self.instancia.vertices.findAt(
            (self.pt_tras_esq,),
            (self.pt_tras_dir,),
            (self.pt_bico_esq,),
            (self.pt_bico_dir,)
        )

        regiao_apoios = self.meu_assembly.Set(
            vertices=vertices_engaste, name='apoios')
        self.meu_modelo.DisplacementBC(name='engaste', createStepName='Initial', region=regiao_apoios,
                                       u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)


def modelar(chassi):
    chassi.geometria()
    chassi.gerar_malha()
    chassi.assembly()
    chassi.carga_engaste()


'''
comp_joelho é definido para o joelho de 25mm, encontrado empiricamente


'''

if __name__ == '__main__':
    chassi_1 = GerarModeloChassi(L=53.5, h=15, b=41, D=2.5, t=0.15, R=11, peso=50, proporcao_banco=0.5, nome = 'Model-1', comp_joelho_cm=3.25)
    modelar(chassi_1)
