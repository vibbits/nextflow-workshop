# QUIZ: Nextflow for reproducible and automated data analysis

1. With what operator can you convert a queue channel to a value channel?
   a. `.view()`
   b. `.first()`
   c. `. collect()`
   d. `.mix()`
   
2. On what language is NextFlow based?
   a. Groovy
   b. Python
   c. Java
   d. R
   e. Perl

3. What is the correct process input notation for this channel entry?

   `[sample1, /data/sample1.fq.gz]`

   a. val(id), path(fq)
   b. tuple val(id), path(fq)
   c. tuple(val(id), path(fq))
   d. val(input)

4. What is the correct process include syntax?
   a. include() to 'main.nf'
   b. include{} from './main.nf'
   c. ./main.nf to {}

5. With what directive can you differentiate  mutiple executions of the same process?

- Open answer in mentimeter or related
  
6.  What is the correct command to download a pipeline?
  a. nextflow run
  b. nextflow config
  c. nextflow pull
  d. nextflow download

7.  With what operator can you manipulate the channel structure
  a. `.map()`
  b. `.view()`
  c. `.mix()`
  d. `.manipulate()`

8. How can you acess ythe value from `--reads`?

   `nextflow run main.nf --reads*.fq.gz`

   a. parameters.reads
   b. reads
   c. reads_param
   d. params.reads

9. with what directive can you differentiate multiple executions of the same process

    - open anser (`tag` should be answerd)

10. What nextflow command do you use to execute a file called 'main.nf'?

      - open answer (`nextflow run main.nf` should be the answer)
12. Which is the better worflow manager?
  a. Nextflow
  b. More Nextflow
  c. How about Nextflow
  d. Nextflow again
  e. Definetely Nextflow
