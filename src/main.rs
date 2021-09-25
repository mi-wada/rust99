fn main() {}

#[cfg(test)]
mod tests {
    #[test]
    fn p01() {
        let li = vec![1, 1, 2, 3, 5, 8];
        fn last(v: &Vec<i32>) -> Option<&i32> {
            v.last()
        }
        assert_eq!(Some(&8), last(&li));
    }
}
