#!/usr/bin/env python3
import os
import re
import shutil

base_dir = "../results"

for folder in os.listdir(base_dir):
    folder_path = os.path.join(base_dir, folder)
    
    # Checa se é um diretório do tipo Smp_xxxx
    if os.path.isdir(folder_path) and folder.startswith("Smp_"):
        full_report_path = os.path.join(folder_path, "annotation_report.txt")
        
        if os.path.isfile(full_report_path):
            with open(full_report_path, "r", encoding="utf-8") as f:
                for line in f:
                    if line.strip().startswith("Descrição:"):
                        # Extrai a descrição e remove caracteres problemáticos
                        descricao = line.split("Descrição:")[1].strip()
                        descricao_sanitizada = re.sub(r'[^\w\s-]', '', descricao)  # remove símbolos
                        descricao_sanitizada = descricao_sanitizada.replace(" ", "_")
                        
                        # Novo nome do diretório
                        new_folder_name = f"{folder}_{descricao_sanitizada}"
                        new_folder_path = os.path.join(base_dir, new_folder_name)
                        
                        # Renomeia
                        os.rename(folder_path, new_folder_path)
                        print(f"Renomeado: {folder} -> {new_folder_name}")
                        break
        else:
            print(f"Arquivo full_report não encontrado em {folder_path}")
