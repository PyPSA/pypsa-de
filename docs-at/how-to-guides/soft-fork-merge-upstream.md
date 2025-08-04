# Soft Fork Upstream Merge 

add the wanted upstream remotes to git
```sh 
git remote add pypsa-de https://github.com/PyPSA/pypsa-de.git
```

fetch the upstream remotes
```sh 
git fetch pypsa-de 
```

checkout a new branch to prevent conflicts in your main branch 
```sh 
git checkout -b pypsa-de-merge
```

merge the upstream branch into your branch 
```sh
git merge pypsa-de/main
```

resolve any merge conflicts in your IDE and assert that the code is working as expected. 
Open a pull request or merge directly into origin/main if you've got permissions to do so. 
```sh
git merge main
```