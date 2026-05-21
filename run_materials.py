# -*- coding: utf-8 -*-
from abaqus import *
from abaqusConstants import *
import job

# --- CONFIGURAÇÕES (Edite aqui) ---
nome_modelo = 'Model-1'


nome_da_secao_solida = 'secao_chassi'  

# Lista: [Nome do Material, Young (Pa), Poisson, Densidade (kg/m3)]
novos_materiais = [
    # Aço (Referência)
    ('Aco',      200e9, 0.30, 7850.0),

    # Alumínio (Mais leve, menos rígido)
    ('Aluminio', 69e9,  0.33, 2700.0),

    # PVC (Muito leve, muito flexível)
    ('PVC',      3e9,   0.38, 1400.0)
]

# Acessa o modelo
myModel = mdb.models[nome_modelo]

# --- INÍCIO DO LOOP ---
for nome_mat, E, Nu, Rho in novos_materiais:

    print('--> Processando: ' + nome_mat)

    # 1. CRIA OU ATUALIZA O MATERIAL
    # Verifica se o material já existe para não dar erro
    if nome_mat not in myModel.materials:
        myMat = myModel.Material(name=nome_mat)
    else:
        myMat = myModel.materials[nome_mat]

    # Define as propriedades (Sobrescreve se já existirem)
    myMat.Density(table=((Rho, ), ))
    myMat.Elastic(table=((E, Nu), ))

    # 2. VINCULA O MATERIAL À SEÇÃO SÓLIDA
    # Aqui é a mágica: Ele pega a sua seção sólida existente e troca o "recheio" dela
    myModel.sections[nome_da_secao_solida].setValues(material=nome_mat)

    # 3. CRIA O JOB
    nome_job = 'Simulacao_' + nome_mat

    # Cria o Job ligado ao modelo (que agora está com o material novo)
    mdb.Job(name=nome_job, model=nome_modelo,
            description='Chassi 3D com '+nome_mat)

    # 4. SUBMETE (RODA) O JOB
    # O script vai mandar rodar e pular para o próximo imediatamente.
    # Se seu PC for fraco, pode ser bom esperar um acabar antes do outro.
    mdb.jobs[nome_job].submit()

    print('    Job ' + nome_job + ' enviado para calculo.')

print('------------------------------------------------')
print('Sucesso! Todos os jobs foram criados e iniciados.')
