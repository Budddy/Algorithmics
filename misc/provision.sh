#!/bin/bash

echo "Going root"
sudo -i

echo "Updating package cache..."
apt-get update > /dev/null

echo "Installing Java..."
apt-get install openjdk-7-jre-headless -y

echo "installing g++"
apt-get install g++ -y

echo "Finished provisioning. Elasticsearch: http://localhost:9200/_plugin/head/; NodeJS: http://localhost:3000"
