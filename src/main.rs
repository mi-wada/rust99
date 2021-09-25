fn main() {}

#[cfg(test)]
mod tests {
    #[test]
    fn p01() {
        let li = vec![1, 1, 2, 3, 5, 8];
        fn last(v: &Vec<i32>) -> Option<i32> {
            v.last().copied()
        }
        assert_eq!(Some(8), last(&li));
    }

    #[test]
    fn p02() {
        let li = vec![1, 1, 2, 3, 5, 8];
        fn penultimate(v: &Vec<i32>) -> Option<i32> {
            v[..v.len() - 1].last().copied()
        }
        assert_eq!(Some(5), penultimate(&li));
    }
}
