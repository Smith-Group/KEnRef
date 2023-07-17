/*
 * KEnRefIinitializer.h
 *
 *      Author: amr
 */

#ifndef KENREFINITIALIZER_H_
#define KENREFINITIALIZER_H_

#include <filesystem>

class KEnRefIinitializer{
public:
	virtual void loadNOEData(std::filesystem::path) = 0;
	virtual void loadStructuralData(std::filesystem::path) = 0;
	virtual ~KEnRefIinitializer() = 0;
};


class DefaultKEnRefInitializer: public KEnRefIinitializer {
public:
	~DefaultKEnRefInitializer() override = 0;
	void loadNOEData(std::filesystem::path) override;
};

#endif /* KENREFINITIALIZER_H_ */
