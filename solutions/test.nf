samples_ch = Channel
                .fromPath('exercises/input.csv')
                .splitCsv(header:true)
                .view{ row -> tuple(row.sampleId, file(row.forward_read), file(row.reverse_read)) }