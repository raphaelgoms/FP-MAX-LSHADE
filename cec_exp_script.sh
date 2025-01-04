#! /bin/bash

while read -r config; do
echo "[$(date)]: Iniciou a exec. para a configuração ${config}" | tee -a log_v2.txt
for func in {1..30}; do
		echo ./solver --cec $func $config 
		./solver --cec $func $config
	done
echo "[$(date)]: Finalizou para a configuração $config" | tee -a log_v2.txt
done < algconfigs