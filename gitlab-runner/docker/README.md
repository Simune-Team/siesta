Steps to setup docker image into gitlab-ci:
```bash
sudo docker build --tag debian-siesta .
sudo docker login
sudo docker tag debian-siesta jalberdi004/debian-siesta
sudo docker push jalberdi004/debian-siesta
```
