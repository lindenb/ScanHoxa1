CLASSPATH=${JVARKIT_DIST}/vcf2table.jar

scanHoxa1.jar: ScanHoxa1.java
	rm -rf TMP
	mkdir TMP
	javac -sourcepath /LAB-DATA/BiRD/users/lindenbaum-p/src/jvarkit/src/main/java:. -d TMP -cp $(CLASSPATH) $<
	jar cvf $@ -C TMP .
	rm -rf TMP
