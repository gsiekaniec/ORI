add_test(
	NAME test_length
	COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/index_creation/; ./test_length.sh"
)

add_test(
	NAME test_makebfQ
	COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/index_creation/; ./test_makebfQ.sh"
)

add_test(
	NAME test_distance
	COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/index_creation/; ./test_distance.sh"
)

add_test(
	NAME test_threshold
	COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/index_creation/; ./test_threshold.sh"
)

add_test(
	NAME test_merge
	COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/index_creation/; ./test_merge.sh"
)

add_test(
	NAME test_clean_merge
	COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/index_creation/; ./test_clean_merge.sh"
)

add_test(
	NAME test_merge_length
	COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/index_creation/; ./test_merge_length.sh"
)

add_test(
	NAME test_cluster
	COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/index_creation/; ./test_cluster.sh"
)

add_test(
	NAME test_build
	COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/index_creation/; ./test_build.sh"
)

add_test(
	NAME test_suppr_bad_reads
	COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/query/; ./test_suppr_bad_reads.sh"
)




