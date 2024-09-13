#!/bin/bash

read -p "Enter version number " VERSION
docker login
echo Building Docker image $VERSION
docker build --platform linux/amd64 -t vkhadka/reader-test:"$VERSION" .
echo Pushing to Docker Hub
docker push vkhadka/reader-test:"$VERSION"